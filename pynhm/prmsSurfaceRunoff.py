from .storageUnit import StorageUnit


class prmsSurfaceRunoff(StorageUnit):
    def __init__(
        self,
        id,
        area,
        forcing,
        verbose,
        hru_percent_imperv,
        imperv_stor_start,
        imperv_stor_max,
        carea_min,
        carea_max,
        soil_rechr_max_frac,
        smidx_coef,
        smidx_exp,
        soil_moisture0_in=None,
    ):
        self.hru_percent_imperv = hru_percent_imperv
        self.imperv_stor_max = imperv_stor_max
        self.imperv_stor_new = imperv_stor_start
        self.imperv_stor_old = None
        self.carea_min = carea_min
        self.carea_max = carea_max
        self.soil_rech_max_frac = (soil_rechr_max_frac,)
        self.smidx_coef = smidx_coef
        self.smidx_exp = smidx_exp
        self.soil_moisture0_in = soil_moisture0_in
        self.soil_moisture0 = 0.0
        super().__init__("sro", id, area, forcing, verbose)
        self.output_column_names = [
            "date",
            "net_precipitation",
            "impervious_runoff",
            "impervious_ds",
            "impervious_et",
            "pervious_runoff",
            "infiltration",
            "residual",
            "impervious_stor_new",
            "impervious_stor_old",
        ]
        return

    def advance(self, itime_step):
        super().advance(itime_step)
        self.imperv_stor_old = self.imperv_stor_new
        if self.soil_moisture0_in is None:
            self.soil_moisture0 = (
                1.0  # todo: this needs to come from the soil zone
            )
        else:
            self.soil_moisture0 = self.soil_moisture0_in[itime_step]
        return

    def calculate(self, time_length, smidx_in=None):

        # Retrieve forcings
        net_precip = self.inflow_volumes["net_precipitation"]
        pot_et = self.forcing.pot_et_current
        snowmelt = 0.0
        upslope_hortonian = 0.0

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
            impervious_runoff = 0.0
            impervious_stor += inflow
        impervious_runoff *= impervious_area

        # storage change in pervious area
        self.imperv_stor_new = impervious_stor
        used_storage = self.imperv_stor_new - self.imperv_stor_old
        used_storage *= impervious_area

        # calculate pervious runoff
        smidx = self.soil_moisture0 + 0.5 * net_precip_length
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
