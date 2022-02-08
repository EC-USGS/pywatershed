"""

"""


class driver:
    def __init__(
        self,
        current_time,
        end_time,
        delta_time,
        forcings,
        storage_units,
        verbose=False,
    ):
        print("Initializing simulation...")
        self.current_time = current_time
        self.end_time = end_time
        self.delta_time = delta_time
        self.forcings = forcings
        self.storage_units = storage_units
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
        self.forcings.advance(self.itime_step, self.current_time)

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
