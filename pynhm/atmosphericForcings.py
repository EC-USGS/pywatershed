class AtmosphericForcings:
    def __init__(self, precip, pot_et):
        self.precip = precip
        self.pot_et = pot_et
        self.pot_et_consumed = None
        self.current_date = None

    def advance(self, itime_step, current_date):
        self.precip_current = self.precip[itime_step]
        self.pot_et_current = self.pot_et[itime_step]
        self.pot_et_consumed = 0.0
        self.current_date = current_date

    def consume_pot_et(self, requested_et):
        et = requested_et
        available_et = self.pot_et_current - self.pot_et_consumed
        if et > available_et:
            et = available_et
        self.pot_et_consumed += et
        return et
