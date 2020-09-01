import Bio.AlignIO
import re
import yaml

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

def text_input(input_file):
    with open(input_file) as f:
        return f.read()

def text_output(string, output_file):
    with open(output_file, 'w') as f:
        f.write(string)

def sequence_input(input_file, format='fasta'):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=format))

    return { record.name: str(record.seq) for record in sequences }

def convert_sequences(input_file, input_format, output_file, output_format):
    with open(input_file) as f:
        sequences = list(Bio.AlignIO.parse(f, format=input_format))

    with open(output_file, 'w') as f:
        Bio.AlignIO.write(sequences, f, output_format)

def yaml_input(input_file):
    with open(input_file) as f:
        return yaml.safe_load(f)

def yaml_output(dict, output_file):
    with open(output_file, 'w') as f:
        yaml.dump(dict, f)