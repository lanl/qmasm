########################################
# Store and filter solutions for QMASM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

class Solutions(object):
    "Represent all near-minimal states of a spin system."

    def __init__(self, answer, problem, all_vars):
        # Store our arguments.
        self.answer = answer
        self.problem = problem
        self.all_vars = all_vars
