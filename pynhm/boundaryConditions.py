class BoundaryConditions:
    def __init__(self):
        self.boundary_conditions = []
        self.current_boundary_conditions = []
        self.current_date = None

    def add_boundary_condition(self, boundary_name, boundary_condition):
        self.boundary_conditions.append((boundary_name, boundary_condition))

    def advance(self, itime_step, current_date):
        self.current_boundary_conditions = []
        for idx, (name, values) in enumerate(self.boundary_conditions):
            self.current_boundary_conditions.append((name, values[itime_step]))
        self.current_date = current_date
