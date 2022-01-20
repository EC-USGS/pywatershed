import os
import numpy as np


class AtmosphericForcings():
    def __init__(self, precip, pot_et):
        self.precip = precip
        self.pot_et = pot_et
        self.pot_et_consumed = None
        self.current_date = None

    def advance(self, itime_step, current_date):
        self.precip_current = self.precip[itime_step]
        self.pot_et_current = self.pot_et[itime_step]
        self.pot_et_consumed = 0
        self.current_date= current_date

    def consume_pot_et(self, requested_et):
        et = requested_et
        available_et = self.pot_et_current - self.pot_et_consumed
        if et > available_et:
            et = available_et
        self.pot_et_consumed += et
        return et


class StorageUnit():

    def __init__(self, storage_type, id, area, forcing, verbose):
        self.storage_type = storage_type
        self.id = id
        self.area = area
        self.forcing = forcing
        self.verbose = verbose
        self.residual_old = None
        self.residual_new = None
        self.inflow_volumes = None
        self.outflow_volumes = None
        self.recipients = []
        self.output_data = []
        self.output_column_names = []
        self.advance()
        return

    def register_recipient(self, recipient, process_name):
        recipient_info = (recipient, process_name)
        self.recipients.append(recipient_info)
        return

    def advance(self):
        self.residual_old = self.residual_new
        self.inflow_volumes = {}
        self.outflow_volumes = {}
        return

    def add_water(self, inflow_name, inflow_rate, time_length):
        self.inflow_volumes[inflow_name] = inflow_rate * time_length
        return

    def remove_water(self, outflow_name, outflow_rate, time_length):
        self.outflow_volumes[outflow_name] = outflow_rate * time_length
        return

    def calculate(self, time_length):
        residual_new = 0.
        for inflow_name in self.inflow_volumes:
            vi = self.inflow_volumes[inflow_name]
            residual_new += vi
        for outflow_name in self.outflow_volumes:
            vo = self.outflow_volumes[outflow_name]
            residual_new += vo
        self.residual_new = residual_new
        return

    def get_budget_summary_str(self):
        s = self.get_name()
        s += f" inflow {self.inflow_volumes} outflow {self.outflow_volumes} " \
             f"residual {self.residual_new} "
        return s

    def get_name(self):
        return f"{self.storage_type}{self.id}"

    def finalize(self):
        return


class Canopy(StorageUnit):

    def __init__(self, id, area, forcing, verbose,
                 intcp_stor_start, intcp_stor_max, covden):
        volume_start = intcp_stor_start * covden
        self.intcp_stor_max = intcp_stor_max
        self.intcp_stor_new = intcp_stor_start
        self.intcp_stor_old = None
        self.covden = covden
        super().__init__("cnp", id, area, forcing, verbose)
        self.output_column_names = ["date", "precip", "aet", "intcp_ds", "net_precip", "residual",
                                    "intcp_stor_new", "intcp_stor_old"]
        return

    def advance(self):
        super().advance()
        self.intcp_stor_old = self.intcp_stor_new
        return

    def calculate(self, time_length):

        # Retrieve forcings
        precip = self.forcing.precip_current
        pot_et = self.forcing.pot_et_current

        # Calculate canopy ET as the minimum of the canopy storage volume
        # and the potential ET.  Then let the forcing module know that
        intcp_stor = self.intcp_stor_old
        requested_et = min(intcp_stor, pot_et)
        canopy_et = self.forcing.consume_pot_et(requested_et)
        intcp_stor -= canopy_et

        # calculate change in storage in canopy and precipitation throughfall
        available_storage = self.intcp_stor_max - intcp_stor
        if precip > available_storage:
            precipitation_throughfall = precip - available_storage
            intcp_stor = self.intcp_stor_max
        else:
            precipitation_throughfall = 0.
            intcp_stor += precip
        self.intcp_stor_new = intcp_stor
        used_storage = (self.intcp_stor_new - self.intcp_stor_old)
        used_storage *= self.area * self.covden

        # net precipitation
        net_precipitation = precip * (1 - self.covden) + precipitation_throughfall * self.covden
        net_precipitation *= self.area

        # add and remove volumetric flow rates
        self.add_water("precipitation", precip * self.area, time_length)
        self.remove_water("aet", -canopy_et * self.area * self.covden, time_length)
        self.remove_water("intcp_ds", -used_storage, time_length)
        self.remove_water("net_precipitation", -net_precipitation, time_length)

        # Send net_precipitation to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "net_precipitation":
                recipient.add_water("net_precipitation", net_precipitation, time_length)

        # Calculate residual and save results
        super().calculate(time_length)
        output = [self.forcing.current_date]
        output.append(precip * self.area)
        output.append(-canopy_et * self.area * self.covden)
        output.append(-used_storage)
        output.append(-net_precipitation)
        output.append(self.residual_new)
        output.append(self.intcp_stor_new)
        output.append(self.intcp_stor_old)
        self.output_data.append(output)
        return


class SurfaceRunoff(StorageUnit):
    def __init__(self, id, area, forcing, verbose, hru_percent_imperv,
                 imperv_stor_start, imperv_stor_max,
                 carea_min, carea_max, soil_rechr_max_frac, smidx_coef, smidx_exp):
        self.hru_percent_imperv = hru_percent_imperv
        self.imperv_stor_max = imperv_stor_max
        self.imperv_stor_new = imperv_stor_start
        self.imperv_stor_old = None
        self.carea_min = carea_min
        self.carea_max = carea_max
        self.soil_rech_max_frac = soil_rechr_max_frac,
        self.smidx_coef = smidx_coef
        self.smidx_exp = smidx_exp
        super().__init__("sro", id, area, forcing, verbose)
        self.output_column_names = ["date", "net_precipitation", "impervious_runoff", "impervious_ds",
                                    "impervious_et", "pervious_runoff", "infiltration", "residual",
                                    "impervious_stor_new", "impervious_stor_old"]
        return

    def advance(self):
        super().advance()
        self.imperv_stor_old = self.imperv_stor_new
        return

    def calculate(self, time_length):

        # Retrieve forcings
        net_precip = self.inflow_volumes["net_precipitation"]
        pot_et = self.forcing.pot_et_current
        snowmelt = 0.
        upslope_hortonian = 0.

        # initialize
        impervious_area = self.area * self.hru_percent_imperv

        # Calculate impervious ET as the minimum of the impervious storage volume
        # and the potential ET.  Then let the forcing module know that
        impervious_stor = self.imperv_stor_old
        requested_et = min(impervious_stor, pot_et)
        impervious_et = self.forcing.consume_pot_et(requested_et)
        impervious_stor -= impervious_et
        impervious_et *= impervious_area

        # calculate impervious runoff
        net_precip_length = net_precip / self.area
        inflow = net_precip_length + snowmelt
        available_storage = self.imperv_stor_max - impervious_stor
        if inflow > available_storage:
            impervious_runoff = inflow - available_storage
            impervious_stor = self.imperv_stor_max
        else:
            impervious_runoff = 0.
            impervious_stor += inflow
        impervious_runoff *= impervious_area

        # storage change in pervious area
        self.imperv_stor_new = impervious_stor
        used_storage = (self.imperv_stor_new - self.imperv_stor_old)
        used_storage *= impervious_area

        # calculate pervious runoff
        smidx = 1. # todo: this needs to come from the soil zone
        ca_fraction = self.smidx_coef * 10 ** (self.smidx_exp * smidx)
        if ca_fraction > self.carea_max:
            ca_fraction = self.carea_max
        pervious_net_precip = net_precip * (1 - self.hru_percent_imperv)
        pervious_runoff = ca_fraction * pervious_net_precip
        infiltration = pervious_net_precip - pervious_runoff

        # add and remove volumetric flow rates
        self.remove_water("impervious_runoff", -impervious_runoff, time_length)
        self.remove_water("impervious_ds", -used_storage, time_length)
        self.remove_water("impervious_et", -impervious_et, time_length)
        self.remove_water("pervious_runoff", -pervious_runoff, time_length)
        self.remove_water("infiltration", -infiltration, time_length)

        # Calculate residual and save results
        super().calculate(time_length)
        output = [self.forcing.current_date]
        output.append(net_precip)
        output.append(-impervious_runoff)
        output.append(-used_storage)
        output.append(-impervious_et)
        output.append(-pervious_runoff)
        output.append(-infiltration)
        output.append(self.residual_new)
        output.append(self.imperv_stor_new)
        output.append(self.imperv_stor_old)

        # append output
        self.output_data.append(output)
        return

class CanopySimple(StorageUnit):

    def __init__(self, id, area, volume_start, volume_max, forcing):
        super().__init__("cnp", id, area, volume_start, forcing)
        self.volume_max = volume_max
        return

    def calculate(self, time_length):

        # add rainfall and update volume stored in canopy
        precip = self.forcing.precip_current
        self.add_water(precip, time_length)
        super().calculate(time_length)

        # Calculate canopy ET as the minimum of the canopy storage volume
        # and the potential ET.  Then let the forcing module know that
        pot_et = self.forcing.pot_et_current
        requested_et = min(self.volume_new, pot_et)
        canopy_et = self.forcing.consume_pot_et(requested_et)
        self.remove_water(-canopy_et, time_length)
        super().calculate(time_length)

        # Calculate net precipitation
        d = self.volume_new - self.volume_max
        net_precipitation = 0.
        if d > 0:
            net_precipitation = d
        self.remove_water(-net_precipitation, time_length)

        # Send net_precipitation to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "net_precipitation":
                recipient.add_water(net_precipitation, time_length)
        super().calculate(time_length)
        return


class ImperviousStorageArea(StorageUnit):
    def __init__(self, id, area, volume_start, volume_max, forcing):
        super().__init__("imp", id, area, volume_start, forcing)
        self.volume_max = volume_max
        return

    def calculate(self, time_length):

        # Calculate impervious ET as the minimum of the impervious storage volume
        # and the potential ET.  Then let the forcing module know that
        pot_et = self.forcing.pot_et_current
        requested_et = min(self.volume_new, pot_et)
        impervious_et = self.forcing.consume_pot_et(requested_et)
        self.remove_water(-impervious_et, time_length)
        super().calculate(time_length)

        # Calculate runoff
        d = self.volume_new - self.volume_max
        runoff = 0.
        if d > 0:
            runoff = d
        self.remove_water(-runoff, time_length)

        # Send runoff to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "runoff":
                recipient.add_water(runoff, time_length)
        super().calculate(time_length)
        return


class PerviousStorageArea(StorageUnit):
    def __init__(self, id, area, volume_start, volume_max, forcing,
                 carea_min, carea_max, activeET=False):
        super().__init__("prv", id, area, volume_start, forcing)
        self.volume_max = volume_max
        self.carea_min = carea_min
        self.carea_max = carea_max
        self.activeET = activeET
        self.soilRechargeZone = None
        return

    def set_soilRechargeZone(self,  soilRechargeZone):
        self.soilRechargeZone = soilRechargeZone
        return

    def calculate(self, time_length):

        # Calculate runoff
        net_precip = 0.
        if len(self.inflow_volumes) > 0:
            # assume net_precip is first entry in inflow_volumes list
            net_precip = self.inflow_volumes[0]

        # Calculate runoff using linear cap equation
        rechr = self.soilRechargeZone.volume_old
        rechr_max = self.soilRechargeZone.volume_max
        ca_fraction = self.carea_min + ((self.carea_max - self.carea_min) * rechr / rechr_max)
        runoff = ca_fraction * net_precip
        self.remove_water(-runoff, time_length)
        super().calculate(time_length)

        # Send runoff to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "runoff":
                recipient.add_water(runoff, time_length)

        # calculate infiltration, which is any remaining water
        infiltration = self.volume_new
        self.remove_water(-infiltration, time_length)

        # Send infiltration to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "infiltration":
                recipient.add_water(infiltration, time_length)

        super().calculate(time_length)
        return


class SoilRechargeZone(StorageUnit):

    def __init__(self, id, area, volume_start, volume_max, forcing, soil2gwmax):
        super().__init__("srz", id, area, volume_start, forcing)
        self.volume_max = volume_max
        self.soil2gwmax = soil2gwmax
        return

    def calculate(self, time_length):

        infiltration = self.inflow_volumes[0]
        soil_rechr = min(self.volume_old + infiltration, self.volume_max)
        self.add_water(soil_rechr, time_length)

        # Calculate soil recharge zone ET as the minimum of the soil recharge zone storage volume
        # and the potential ET.  Then let the forcing module know that
        pot_et = self.forcing.pot_et_current
        requested_et = min(self.volume_new, pot_et)
        actual_et = self.forcing.consume_pot_et(requested_et)
        self.remove_water(-actual_et, time_length)
        super().calculate(time_length)

        # Send excess infiltration to any registered recipients
        excess_infiltration = infiltration - soil_rechr
        soil2gwr = min(excess_infiltration, self.soil2gwmax)
        soil2grav = excess_infiltration - soil2gwr

        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "soil2gwr":
                recipient.add_water(soil2gwr, time_length)
            elif process_name == "soil2grav":
                recipient.add_water(soil2grav, time_length)
        return


if __name__ == "__main__":

    print("standalone run of {}".format(os.path.basename(__file__)))

    # initialize
    print("Initializing...")
    rainfall = [4.5, 0., 0., 0., 0.]
    potential_evapotranspiration = len(rainfall) * [1.]

    forcings = AtmosphericForcings(rainfall, potential_evapotranspiration)
    cnp = Canopy(id=0, area=1., volume_start=0., volume_max=2.0, forcing=forcings)
    imp = ImperviousStorageArea(id=1, area=1., volume_start=0., volume_max=0.1, forcing=forcings)
    prv = PerviousStorageArea(id=2, area=1., volume_start=0., volume_max=0.5, forcing=forcings,
                              carea_max=0.75, carea_min=0.25, activeET=False)
    srz = SoilRechargeZone(id=3, area=1., volume_start=0., volume_max=0.5, forcing=forcings, soil2gwmax=1.)
    #gwr = GroundwaterReservoir(id=4, ...)
    #grv = GravityReservoir(id=5, ...)

    # The pervious storage area needs access to rech and rechmax to calculate infiltration
    prv.set_soilRechargeZone(srz)


    # register recipients
    cnp.register_recipient(imp, "net_precipitation")
    cnp.register_recipient(prv, "net_precipitation")
    prv.register_recipient(srz, "infiltration")

    #srz.register_recipient(gwr, "soil2gwr")
    #srz.register_recipient(grv, "soil2grav")

    # The order of this list of storage units should ultimately be based on a dag
    storage_units = [cnp, imp, prv, srz]

    # do timesteps
    time_length = 1.0

    for itime_step in range(len(rainfall)):

        # set precip and pot_et
        precip = rainfall[itime_step]
        pot_et = potential_evapotranspiration[itime_step]
        print(f"\nProcessing step {itime_step} with rainfall = {precip} and pot_et = {pot_et}")

        forcings.advance(itime_step)

        for storage_unit in storage_units:
            storage_unit.advance()

        for storage_unit in storage_units:
            storage_unit.calculate(time_length)

        for storage_unit in storage_units:
            print(storage_unit.get_budget_summary_str())

    print("Finalizing...")
    for storage_unit in storage_units:
        storage_unit.finalize()
