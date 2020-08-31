import Bio.AlignIO
import re

def update_dict(dict, **kwargs):
    res = dict.copy()
    res.update(**kwargs)
    return res

def cmd_kwargs(kwarg_prefix="-", **kwargs):
    pairs = [(kwarg_prefix + key, str(value)) for key, value in kwargs.items()]
    return [x for pair in pairs for x in pair]

def build_shell_command(path, args, kwargs, arg_prefix="--", kwarg_prefix="-"):
    return " ".join(
        [path] +
        [arg_prefix + str(arg) for arg in args] +
        cmd_kwargs(kwarg_prefix, **kwargs)
    )

def convert_sequences(input_file, input_format, output_file, output_format):
    with open(input_file) as f:
        sequences = list(Bio.AlignIO.parse(f, format=input_format))

    with open(output_file, 'w') as f:
        Bio.AlignIO.write(sequences, f, output_format)