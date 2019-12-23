###################################
# Output QUBOs in various formats #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import datetime
import json
import random
import sys

class OutputMixin(object):
    "Provide functions for outputting problems and solutions."

    def open_output_file(self, oname):
        "Open a file or standard output."
        if oname == "<stdout>":
            outfile = sys.stdout
        else:
            try:
                outfile = open(oname, "w")
            except IOError:
                self.abend('Failed to open %s for output' % oname)
        return outfile

    def output_qmasm(self, outfile):
        "Output weights and strengths as a flattened QMASM source file."
        for p in self.program:
            outfile.write("%s\n" % p.as_str())

    def output_bqpjson(self, outfile, as_qubo, problem):
        "Output weights and strengths in bqpjson format, either Ising or QUBO."
        # Prepare the "easy" fields.
        bqp = {}
        bqp["version"] = "1.0.0"
        bqp["id"] = random.randint(2**20, 2**60)
        bqp["scale"] = 1.0
        bqp["offset"] = 0.0
        if as_qubo:
            bqp["variable_domain"] = "boolean"
        else:
            bqp["variable_domain"] = "spin"

        # Prepare the list of all variables.
        var_ids = set(problem.weights.keys())
        for q1, q2 in problem.strengths.keys():
            var_ids.add(q1)
            var_ids.add(q2)
        bqp["variable_ids"] = sorted(var_ids)

        # Prepare the linear terms.
        lin_terms = []
        for q, wt in sorted(problem.weights.items()):
            lin_terms.append({
                "id": q,
                "coeff": wt})
        bqp["linear_terms"] = lin_terms

        # Prepare the quadratic terms.
        quad_terms = []
        strengths = self.canonicalize_strengths(problem.strengths)
        for (q1, q2), wt in sorted(strengths.items()):
            quad_terms.append({
                "id_tail": q1,
                "id_head": q2,
                "coeff": wt})
        bqp["quadratic_terms"] = quad_terms

        # Prepare some metadata.
        metadata = {}
        if as_qubo:
            metadata["description"] = "QUBO problem compiled by QMASM (https://github.com/lanl/qmasm)"
        else:
            metadata["description"] = "Ising problem compiled by QMASM (https://github.com/lanl/qmasm)"
        metadata["command_line"] = self.get_command_line()
        metadata["generated"] = datetime.datetime.utcnow().isoformat()
        if hasattr(problem, "embedding"):
            # Physical problem
            def attempt_assign(key, func):
                "Try assigning a key, but don't complain if we can't."
                try:
                    metadata[key] = func()
                except KeyError:
                    pass
            attempt_assign("dw_url", lambda: os.environ["DW_INTERNAL__HTTPLINK"])
            attempt_assign("dw_solver_name", lambda: self.solver_name)
            props = self.solver.properties
            attempt_assign("dw_chip_id", lambda: props["chip_id"])
            L, M, N = self.chimera_topology(self.solver)
            metadata["chimera_cell_size"] = L*2
            metadata["chimera_degree"] = max(M, N)
            metadata["equivalent_ids"] = sorted(problem.chains)
            metadata["variable_names"] = {s: problem.embedding[n]
                                          for s, n in self.sym_map.symbol_number_items()}
        else:
            metadata["variable_names"] = {s: [n]
                                          for s, n in self.sym_map.symbol_number_items()}
        bqp["metadata"] = metadata

        # Output the problem in JSON format.
        outfile.write(json.dumps(bqp, indent=2, sort_keys=True) + "\n")

    def write_output(self, problem, oname, oformat, as_qubo):
        "Write an output file in one of a variety of formats."

        # Open the output file.
        outfile = self.open_output_file(oname)

        # Output the weights and strengths in the specified format.
        if oformat == "qubist":
            self.output_qubist(outfile, as_qubo, problem)
        elif oformat == "dw":
            self.output_dw(outfile, problem)
        elif oformat == "qbsolv":
            self.output_qbsolv(outfile, problem)
        elif oformat == "qmasm":
            self.output_qmasm(outfile)
        elif oformat == "minizinc":
            self.output_minizinc(outfile, problem)
        elif oformat == "bqpjson":
            self.output_bqpjson(outfile, as_qubo, problem)

        # Close the output file.
        if oname != "<stdout>":
            outfile.close()
