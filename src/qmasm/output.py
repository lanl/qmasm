###################################
# Output QUBOs in various formats #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import datetime
import json
import random
import shlex
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

    def output_minizinc(self, outfile, as_qubo, problem, energy=None):
        "Output weights and strengths as a MiniZinc constraint problem."
        # Write some header information.
        outfile.write("""% Use MiniZinc to minimize a given Hamiltonian.
    %
    % Producer:     QMASM (https://github.com/lanl/qmasm/)
    % Author:       Scott Pakin (pakin@lanl.gov)
    """)
        outfile.write("%% Command line: %s\n\n" % " ".join([shlex.quote(a) for a in sys.argv]))

        # The model is easier to express as a QUBO so convert to that format.
        if as_qubo:
            qprob = problem
        else:
            qprob = problem.convert_to_qubo()

        # Map each qubit to one or more symbols.
        num2syms = {}
        for s, n in self.sym_map.symbol_number_items():
            try:
                # Physical problem
                for pn in qprob.embedding[n]:
                    try:
                        num2syms[pn].append(s)
                    except KeyError:
                        num2syms[pn] = [s]
            except AttributeError:
                # Logical problem
                try:
                    num2syms[n].append(s)
                except KeyError:
                    num2syms[n] = [s]
        for n in num2syms.keys():
            num2syms[n].sort(key=lambda s: ("$" in s, s))

        # Find the character width of the longest list of symbol names.
        max_sym_name_len = max([len(repr(ss)) - 1 for ss in num2syms.values()] + [7])

        # Output all QMASM variables as MiniZinc variables.
        qubits_used = set(qprob.weights.keys())
        qubits_used.update([qs[0] for qs in qprob.strengths.keys()])
        qubits_used.update([qs[1] for qs in qprob.strengths.keys()])
        for q in sorted(qubits_used):
            outfile.write("var 0..1: q%d;  %% %s\n" % (q, " ".join(num2syms[q])))
        outfile.write("\n")

        # Define variables representing products of QMASM variables.  Constrain
        # the product variables to be the products.
        outfile.write("% Define p_X_Y variables and constrain them to be the product of qX and qY.\n")
        for q0, q1 in sorted(qprob.strengths.keys()):
            pstr = "p_%d_%d" % (q0, q1)
            outfile.write("var 0..1: %s;\n" % pstr)
            outfile.write("constraint %s >= q%d + q%d - 1;\n" % (pstr, q0, q1))
            outfile.write("constraint %s <= q%d;\n" % (pstr, q0))
            outfile.write("constraint %s <= q%d;\n" % (pstr, q1))
        outfile.write("\n")

        # Express energy as one, big Hamiltonian.
        scale_to_int = lambda f: int(round(10000.0*f))
        outfile.write("var int: energy =\n")
        weight_terms = ["%8d * q%d" % (scale_to_int(wt), q) for q, wt in sorted(qprob.weights.items())]
        strength_terms = ["%8d * p_%d_%d" % (scale_to_int(s), qs[0], qs[1]) for qs, s in sorted(qprob.strengths.items())]
        all_terms = weight_terms + strength_terms
        outfile.write("  %s;\n" % " +\n  ".join(all_terms))

        # Because we can't both minimize and enumerate all solutions, we
        # normally do only the former with instructions for the user on how to
        # switch to the latter.  However, if an energy was specified, comment
        # out the minimization step and uncomment the enumeration step.
        outfile.write("\n")
        outfile.write("% First pass: Compute the minimum energy.\n")
        if energy == None:
            outfile.write("solve minimize energy;\n")
        else:
            outfile.write("% solve minimize energy;\n")
        outfile.write("""
%% Second pass: Find all minimum-energy solutions.
%%
%% Once you've solved for minimum energy, comment out the "solve minimize
%% energy" line, plug the minimal energy value into the following line,
%% uncomment it and the "solve satisfy" line, and re-run MiniZinc, requesting
%% all solutions this time.  The catch is that you need to use the raw
%% energy value so be sure to modify the output block to show(energy)
%% instead of show(energy/%.10g + %.10g).
""" % (self.minizinc_scale_factor, qprob.bqm.offset))
        if energy == None:
            outfile.write("%constraint energy = -12345;\n")
            outfile.write("%solve satisfy;\n\n")
        else:
            outfile.write("constraint energy = %d;\n" % energy)
            outfile.write("solve satisfy;\n\n")

        # Output code to show the results symbolically.  We output in the same
        # format as QMASM normally does.  Unfortunately, I don't know how to get
        # MiniZinc to output the current solution number explicitly so I had to
        # hard-wire "Solution #1".
        outfile.write("output [\n")
        outfile.write('  "Solution #1 (energy = ", show(energy/%.10g + %.10g), ", tally = 1)\\n\\n",\n' % (self.minizinc_scale_factor, qprob.bqm.offset))
        outfile.write('  "    %-*s  Spin  Boolean\\n",\n' % (max_sym_name_len, "Name(s)"))
        outfile.write('  "    %s  ----  -------\\n",\n' % ("-" * max_sym_name_len))
        outlist = []
        for n, ss in num2syms.items():
            if ss == []:
                continue
            syms = " ".join(ss)
            line = ""
            line += '"    %-*s  ", ' % (max_sym_name_len, syms)
            if as_qubo:
                line += 'show_int(4, q%d), ' % n
            else:
                line += 'show_int(4, 2*q%d - 1), ' % n
            line += '"  ", if show(q%d) == "1" then "True" else "False" endif, ' % n
            line += '"\\n"'
            outlist.append(line)
        outlist.sort()
        outfile.write("  %s\n];\n" % ",\n  ".join(outlist))

    def output_qbsolv(self, outfile, problem):
        "Output weights and strengths in qbsolv format."
        # Determine the list of nonzero weights and strengths.
        qprob = problem.convert_to_qubo()
        output_weights, output_strengths = qprob.weights, qprob.strengths
        max_node = max(list(output_weights.keys()) + [max(qs) for qs in output_strengths.keys()])
        num_nonzero_weights = len([q for q, wt in output_weights.items() if wt != 0.0])
        num_nonzero_strengths = len([qs for qs, wt in output_strengths.items() if wt != 0.0])

        # Assign dummy qubit numbers to qubits whose value is known a priori.
        try:
            n_known = len(qprob.known_values)
        except TypeError:
            n_known = 0
        try:
            extra_nodes = dict(zip(sorted(qprob.known_values.keys()),
                                   range(max_node + 1, max_node + 1 + n_known)))
        except AttributeError:
            extra_nodes = {}
        max_node += n_known
        num_nonzero_weights += n_known
        output_weights.update({num: qprob.known_values[sym]*qmasm.pin_weight
                               for sym, num in extra_nodes.items()})
        sym2num = dict(self.sym_map.symbol_number_items())
        sym2num.update(extra_nodes)

        # Output a name-to-number map as header comments.
        key_width = 0
        val_width = 0
        items = []
        for s, n in sym2num.items():
            if len(s) > key_width:
                key_width = len(s)

            # Map logical to physical if possible.
            try:
                # Physical problem
                known_values = qprob.logical.merged_known_values()
                pin_map = {k: v for k, v in qprob.logical.pinned}
            except AttributeError:
                # Logical problem
                known_values = qprob.merged_known_values()
                pin_map = {k: v for k, v in qprob.pinned}
            try:
                nstr = " ".join([str(n) for n in sorted(qprob.embedding[n])])
            except AttributeError:
                # Logical problem
                nstr = str(n)
            except KeyError:
                try:
                    nstr = "[Pinned to %s]" % repr(pin_map[n])
                except KeyError:
                    try:
                        nstr = "[Provably %s]" % known_values[n]
                    except KeyError:
                        try:
                            same = qprob.logical.contractions[n]
                            nstr = "[Same as %s]" % " ".join(self.sym_map.to_symbols(same))
                        except KeyError:
                            nstr = "[Disconnected]"
            if len(nstr) > val_width:
                val_width = len(nstr)
            items.append((s, nstr))
        items.sort()
        for s, nstr in items:
            outfile.write("c %-*s --> %-*s\n" % (key_width, s, val_width, nstr))

        # Output all nonzero weights and strengths.
        outfile.write("p qubo 0 %d %d %d\n" % (max_node + 1, num_nonzero_weights, num_nonzero_strengths))
        for q, wt in sorted(output_weights.items()):
            if wt != 0.0:
                outfile.write("%d %d %.10g\n" % (q, q, wt))
        for qs, wt in sorted(output_strengths.items()):
            if wt != 0.0:
                outfile.write("%d %d %.10g\n" % (qs[0], qs[1], wt))

    def output_qubist(self, outfile, as_qubo, problem, sampler):
        "Output weights and strengths in Qubist format, either Ising or QUBO."
        # Convert the problem to Ising, scale it for the hardware, then convert
        # to QUBO if requested.
        prob = problem.convert_to_ising()
        prob.autoscale_coefficients(sampler)
        if as_qubo:
            prob = prob.convert_to_qubo()
        output_weights = prob.weights
        output_strengths = prob.strengths

        # Format all weights and all strengths in Qubist format.
        data = []
        for q, wt in sorted(output_weights.items()):
            if wt != 0.0:
                data.append("%d %d %.10g" % (q, q, wt))
        for sp, str in sorted(output_strengths.items()):
            if str != 0.0:
                sp = sorted(sp)
                data.append("%d %d %.10g" % (sp[0], sp[1], str))

        # Output the header and data in Qubist format.
        try:
            num_qubits = sampler.sampler.properties["num_qubits"]
        except KeyError:
            # If the solver lacks a fixed hardware representation we assert
            # that the number of hardware qubits is exactly the number of
            # qubits we require.
            num_qubits = len(output_weights)
        outfile.write("%d %d\n" % (num_qubits, len(data)))
        for d in data:
            outfile.write("%s\n" % d)

    def output_ocean(self, outfile, as_qubo, problem, sampler, sampler_args):
        "Output weights and strengths as an Ocean program, either Ising or QUBO."
        # Select each variable's alphabetically first symbol, favoring symbols
        # without dollar signs.
        all_nums = self.sym_map.all_numbers()
        num2sym = {}
        for n in all_nums:
            syms = list(self.sym_map.to_symbols(n))
            syms.sort(key=lambda s: ("$" in s, s))
            num2sym[n] = syms[0]

        # Output some Python boilerplate.
        outfile.write("#! /usr/bin/env python\n\n")
        outfile.write("import dimod\n")
        physical = getattr(problem, "logical", None) != None
        if physical:
            outfile.write("from dwave.system import DWaveSampler\n\n")
        else:
            outfile.write("from dwave.system import DWaveSampler, EmbeddingComposite\n\n")

        # Output code to set up the problem.
        if physical:
            linear = ", ".join(["%d: %.16g" % e for e in sorted(problem.bqm.linear.items())])
            quadratic = ", ".join(["(%d, %d): %.20g" % (ns[0], ns[1], wt) for ns, wt in sorted(problem.bqm.quadratic.items())])
        else:
            linear = ", ".join(sorted(["'%s': %.20g" % (num2sym[n], wt) for n, wt in problem.bqm.linear.items()]))
            quadratic = ", ".join(sorted(["('%s', '%s'): %.20g" % (num2sym[ns[0]], num2sym[ns[1]], wt) for ns, wt in problem.bqm.quadratic.items()]))
        outfile.write("linear = {%s}\n" % linear)
        outfile.write("quadratic = {%s}\n" % quadratic)
        if as_qubo:
            vtype = "BINARY"
        else:
            vtype = "SPIN"
        outfile.write("bqm = dimod.BinaryQuadraticModel(linear, quadratic, %.5g, dimod.%s)\n\n" % (problem.bqm.offset, vtype))

        # Modify the sampler arguments as necessary.
        sampler_args = {k: v for k, v in sampler_args.items() if v != None}
        try:
            if sampler_args["anneal_schedule"] != None:
                del sampler_args["annealing_time"]
        except KeyError:
            pass
        try:
            pp = sampler_args["postprocess"]
            if pp == "sample":
                sampler_args["postprocess"] = "sampling"
            elif pp == "opt":
                sampler_args["postprocess"] = "optimization"
            else:
                del sampler_args["postprocess"]
        except KeyError:
            pass
        try:
            if sampler_args["num_spin_reversal_transforms"] == 0:
                del sampler_args["num_spin_reversal_transforms"]
        except:
            pass
        arg_str = ", ".join(["%s=%s" % (k, repr(v)) for k, v in sampler_args.items()])

        # Output code to solve the problem.
        if physical:
            outfile.write("sampler = DWaveSampler()\n")
        else:
            outfile.write("sampler = EmbeddingComposite(DWaveSampler())\n")
        outfile.write("result = sampler.sample(bqm, %s)\n" % arg_str)

        # Output code to display the results QMASM-style.
        outfile.write(r'''
data = result.data(fields=["sample", "energy", "num_occurrences"])
wd = max([8] + [len(v) for v in result.variables])
vnames = sorted(result.variables, key=lambda v: ("$" in v, v))
for i in range(len(result.samples())):
    if i > 0:
        print("")
    s = next(data)
    print("Solution #%d (energy = %.4g, tally = %d):\n" % (i + 1, s.energy, s.num_occurrences))
    print("    %-*s  Value" % (wd, "Variable"))
    print("    %s  -----" % ("-"*wd))
    for v in vnames:
        print("    %-*s  %s" % (wd, v, s.sample[v] == 1))
''')

    def write_output(self, problem, oname, oformat, as_qubo, sampler, sampler_args):
        "Write an output file in one of a variety of formats."

        # Open the output file.
        outfile = self.open_output_file(oname)

        # Output the weights and strengths in the specified format.
        if oformat == "qubist":
            self.output_qubist(outfile, as_qubo, problem, sampler)
        if oformat == "ocean":
            self.output_ocean(outfile, as_qubo, problem, sampler, sampler_args)
        elif oformat == "qbsolv":
            self.output_qbsolv(outfile, problem)
        elif oformat == "qmasm":
            self.output_qmasm(outfile)
        elif oformat == "minizinc":
            self.output_minizinc(outfile, as_qubo, problem)
        elif oformat == "bqpjson":
            self.output_bqpjson(outfile, as_qubo, problem)

        # Close the output file.
        if oname != "<stdout>":
            outfile.close()
