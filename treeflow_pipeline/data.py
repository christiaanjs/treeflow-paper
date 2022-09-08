from Bio.SeqRecord import SeqRecord
import Bio.AlignIO
import datetime


def reformat_record_date(record: SeqRecord, date_format="%Y/%m/%d", min_year=1900):
    split = record.id.split("_")
    date = datetime.datetime.strptime(split[-1], date_format).date()
    date_decimal = (float(date.strftime("%j")) - 1) / 366 + float(date.strftime("%Y"))
    new_date_string = f"{date_decimal:.3f}"
    new_id = "_".join(split[:-1] + [new_date_string])
    return SeqRecord(
        seq=record.seq, id=new_id, name=record.name, description=record.description
    )


def convert_dates_to_numeric(
    input_file, input_format, output_file, output_format, date_format="%Y/%m/%d"
):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    new_msa = Bio.AlignIO.MultipleSeqAlignment(
        [reformat_record_date(record, date_format=date_format) for record in sequences]
    )

    Bio.AlignIO.write([new_msa], output_file, output_format)
