from .core_parameter import CoreParameter


class EnsoDiagsParameter(CoreParameter):
    def __init__(self):
        super(EnsoDiagsParameter, self).__init__()
        self.granulate.remove('seasons')
        self.nino_region = 'NINO34'
        self.print_statements = False
        # A list of the reference names to run the diags on.
        self.ref_names = []
        self.ref_timeseries_input = True
        self.test_timeseries_input = True


    def check_values(self):
        if not self.ref_names:
            msg = 'You have no value for ref_names. Caculate test data only'
            print(msg)

        valid_nino_regions = ['NINO3', 'NINO34', 'NINO4']
        if self.nino_region not in valid_nino_regions:
            msg = 'nino_region={} not in {}'.format(
                self.nino_region, valid_nino_regions) 
            raise RuntimeError(msg)

        if not (hasattr(self, 'start_yr') and hasattr(self, 'end_yr')):
            msg = "You need to define both the 'start_yr' and 'end_yr' parameter."
            raise RuntimeError(msg)
