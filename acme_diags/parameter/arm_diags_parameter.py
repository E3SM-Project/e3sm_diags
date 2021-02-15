from .core_parameter import CoreParameter


class ARMDiagsParameter(CoreParameter):
    def __init__(self):
        super().__init__()
        # A list of the reference names to run the diags on.
        self.granulate.remove('seasons')

#    def check_values(self):
#        if not self.ref_names:
#            msg = 'You must have a value for ref_names.'
#            raise RuntimeError(msg)
