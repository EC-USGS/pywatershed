from ..base.storageUnit import StorageUnit


class prmsCanopy(StorageUnit):
    def __init__(
        self,
        id,
        area,
        forcing,
        verbose,
        intcp_stor_start,
        intcp_stor_max,
        covden,
    ):
        volume_start = intcp_stor_start * covden
        self.intcp_stor_max = intcp_stor_max
        self.intcp_stor_new = intcp_stor_start
        self.intcp_stor_old = None
        self.covden = covden
        super().__init__("cnp", id, area, forcing, verbose)
        self.output_column_names = [
            "date",
            "precip",
            "aet",
            "intcp_ds",
            "net_precip",
            "residual",
            "intcp_stor_new",
            "intcp_stor_old",
        ]
        return

    def advance(self, itime_step):
        super().advance(itime_step)
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
            precipitation_throughfall = 0.0
            intcp_stor += precip
        self.intcp_stor_new = intcp_stor
        used_storage = self.intcp_stor_new - self.intcp_stor_old
        used_storage *= self.area * self.covden

        # net precipitation
        net_precipitation = (
            precip * (1 - self.covden)
            + precipitation_throughfall * self.covden
        )
        net_precipitation *= self.area

        # add and remove volumetric flow rates
        self.add_water("precipitation", precip * self.area, time_length)
        self.remove_water(
            "aet", -canopy_et * self.area * self.covden, time_length
        )
        self.remove_water("intcp_ds", -used_storage, time_length)
        self.remove_water("net_precipitation", -net_precipitation, time_length)

        # Send net_precipitation to any registered recipients
        for recipient_info in self.recipients:
            recipient, process_name = recipient_info
            if process_name == "net_precipitation":
                recipient.add_water(
                    "net_precipitation", net_precipitation, time_length
                )

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
