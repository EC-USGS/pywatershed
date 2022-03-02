"""

"""
import datetime

from .forcings.atmosphericForcings import AtmosphericForcings
from .base.boundaryForcings import BoundaryForcings


class driver:
    def __init__(
        self,
        current_time: datetime,
        end_time: datetime,
        delta_time: datetime.timedelta,
        storage_units: list,
        climate_forcings: AtmosphericForcings = None,
        external_forcings: BoundaryForcings = None,
        verbose: int = False,
    ):
        print("Initializing simulation...")
        self.current_time = current_time
        self.end_time = end_time
        self.delta_time = delta_time
        self.climate_forcings = climate_forcings
        self.storage_units = storage_units
        self.external_forcings = external_forcings
        self.verbose = verbose

        self.time_length = float(delta_time.days)
        self.itime_step = 0
        return

    def run(self):
        """Run the simulation"""
        while self.current_time <= self.end_time:

            if self.verbose:
                print(f"processing...{self.current_time}")

            # run the current time
            self.update()

        # finalize the simulation
        self.finalize()

        return

    def update(self):
        """Run a time step"""
        if self.climate_forcings is not None:
            self.climate_forcings.advance(self.itime_step, self.current_time)
        if self.external_forcings is not None:
            self.external_forcings.advance(self.itime_step, self.current_time)

        for storage_unit in self.storage_units:
            storage_unit.advance(self.itime_step)

        for storage_unit in self.storage_units:
            storage_unit.calculate(self.time_length)

        if self.verbose:
            for storage_unit in self.storage_units:
                print(storage_unit.get_budget_summary_str())

        self.current_time += self.delta_time
        self.itime_step += 1

    def finalize(self):
        print("Finalizing simulation...")
        for storage_unit in self.storage_units:
            storage_unit.finalize()
