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
        metadata["variable_names"] = {s: [n]
                                      for s, n in self.sym_map.symbol_number_items()}
        bqp["metadata"] = metadata

        # Output the problem in JSON format.
        outfile.write(json.dumps(bqp, indent=2, sort_keys=True) + "\n")

    def output_qubist(self, outfile, as_qubo, problem, sampler):
        "Output weights and strengths in Qubist format, either Ising or QUBO."
        if as_qubo and not problem.qubo:
            qprob = problem.convert_to_qubo()
            output_weights, output_strengths = qprob.weights, qprob.strengths
        elif not as_qubo and problem.qubo:
            iprob = problem.convert_to_ising()
            output_weights, output_strengths = iprob.weights, iprob.strengths
        else:
            output_weights = problem.weights
            output_strengths = problem.strengths
        data = []
        for q, wt in sorted(output_weights.items()):
            if wt != 0.0:
                data.append("%d %d %.10g" % (q, q, wt))
        for sp, str in sorted(output_strengths.items()):
            if str != 0.0:
                data.append("%d %d %.10g" % (sp[0], sp[1], str))

        # Output the header and data in Qubist format.
        try:
            num_qubits = sampler.sampler.properties["num_qubits"]
        except KeyError:
            # The Ising heuristic solver is an example of a solver that lacks a
            # fixed hardware representation.  We therefore assert that the number
            # of qubits is exactly the number of qubits we require.
            num_qubits = len(output_weights)
        outfile.write("%d %d\n" % (num_qubits, len(data)))
        for d in data:
            outfile.write("%s\n" % d)

    def write_output(self, problem, oname, oformat, as_qubo, sampler):
        "Write an output file in one of a variety of formats."

        # Open the output file.
        outfile = self.open_output_file(oname)

        # Output the weights and strengths in the specified format.
        if oformat == "qubist":
            self.output_qubist(outfile, as_qubo, problem, sampler)
        elif oformat == "qmasm":
            self.output_qmasm(outfile)
        elif oformat == "minizinc":
            self.output_minizinc(outfile, problem)
        elif oformat == "bqpjson":
            self.output_bqpjson(outfile, as_qubo, problem)

        # Close the output file.
        if oname != "<stdout>":
            outfile.close()
