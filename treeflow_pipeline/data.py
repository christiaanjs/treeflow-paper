from multiprocessing.sharedctypes import Value
import typing as tp
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.AlignIO
import Bio.SeqIO
import datetime
import xml.etree.ElementTree


class RecordDateError(ValueError):
    pass


def reformat_record_date(
    record: SeqRecord,
    date_format="%Y/%m/%d",
    split_char="_",
    split_index=-1,
    date_join_char="_",
    replace_old_date=True,
):
    split = record.id.split(split_char)
    try:
        date = datetime.datetime.strptime(split[split_index], date_format).date()
    except ValueError as ex:
        raise RecordDateError(ex)
    date_decimal = (float(date.strftime("%j")) - 1) / 366 + float(date.strftime("%Y"))
    new_date_string = f"{date_decimal:.3f}"
    if replace_old_date:
        old_components = split[:split_index] + split[split_index + 1 :]
    else:
        old_components = split
    new_id = date_join_char.join([split_char.join(old_components), new_date_string])
    return SeqRecord(seq=record.seq, id=new_id, description="")


def convert_dates_to_numeric(
    input_file,
    input_format,
    output_file,
    output_format,
    date_format="%Y/%m/%d",
    split_char="_",
    split_index=-1,
    date_join_char="_",
    replace_old_date=True,
    exclude_bad_dates=False,
):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    output_records = []
    for record in sequences:
        try:
            output_records.append(
                reformat_record_date(
                    record,
                    date_format=date_format,
                    split_char=split_char,
                    split_index=split_index,
                    date_join_char=date_join_char,
                    replace_old_date=replace_old_date,
                )
            )
        except RecordDateError as ex:
            if not exclude_bad_dates:
                raise ex
            else:
                print(f"Excluding record {record.id} for bad date")
    new_msa = Bio.AlignIO.MultipleSeqAlignment(output_records)

    Bio.AlignIO.write([new_msa], output_file, output_format)


def parse_sequence_value(tag):
    if "value" in tag.attrib:
        return tag.attrib["value"]
    else:
        text = tag.text.strip()
        children = list(tag)
        i = 0
        n = len(children)
        while not text and i < n:
            text = children[i].tail
            i += 1
        return "".join(text.split())


def parse_taxon_name(
    tag, remove_spaces=False, date_trait_dict: tp.Optional[tp.Dict[str, float]] = None
):
    if "taxon" in tag.attrib:
        raw_taxon_name: str = tag.attrib["taxon"]
    else:
        taxon_tag = tag.find("./taxon")
        if "idref" in taxon_tag.attrib:
            raw_taxon_name: str = taxon_tag.attrib["idref"]
        else:
            raise ValueError("Can't find taxon in sequence tag")
    if remove_spaces:
        formatted_name = raw_taxon_name.replace(" ", "_")
    else:
        formatted_name = raw_taxon_name
    if date_trait_dict is None:
        return formatted_name
    else:
        return f"{formatted_name}_{date_trait_dict[raw_taxon_name]}"


DEFAULT_DATE_TRAIT_NAME = "date-forward"


def extract_date_trait_dict(
    xml_root: xml.etree.ElementTree.ElementTree, date_trait_name=DEFAULT_DATE_TRAIT_NAME
) -> tp.Dict[str, float]:
    trait_set_elements = xml_root.findall(f"//*[@traitname='{date_trait_name}']")
    assert (
        len(trait_set_elements) == 1
    ), f"Expected one date trait set, found {len(trait_set_elements)}"
    date_trait_string = trait_set_elements[0].attrib["value"]
    name_date_pairs = [element.split("=") for element in date_trait_string.split(",")]
    return {name.strip(): float(date_string) for name, date_string in name_date_pairs}


def extract_xml_sequences(
    input_file,
    output_file,
    output_format,
    reformat_taxon_name=False,
    date_from_trait_string=False,
    date_trait_name=DEFAULT_DATE_TRAIT_NAME,
):
    seq_xml_root = xml.etree.ElementTree.parse(input_file)
    if date_from_trait_string:
        date_trait_dict = extract_date_trait_dict(
            seq_xml_root, date_trait_name=date_trait_name
        )
    else:
        date_trait_dict = None
    records = [
        SeqRecord(
            Seq(parse_sequence_value(tag)),
            parse_taxon_name(
                tag, remove_spaces=reformat_taxon_name, date_trait_dict=date_trait_dict
            ),
            description="",
        )
        for tag in seq_xml_root.findall(".//sequence")
    ]
    with open(output_file, "w") as f:
        Bio.SeqIO.write(records, f, output_format)


def remove_identical_sequences(input_file, input_format, output_file, output_format):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))

    print(f"Input: {len(sequences)} sequences")

    sequence_dict = {record.seq: record for record in sequences}
    reduced_sequences = list(sequence_dict.values())

    print(f"Output: {len(reduced_sequences)} sequences")
    with open(output_file, "w") as f:
        Bio.SeqIO.write(reduced_sequences, f, output_format)
