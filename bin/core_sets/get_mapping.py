import sys
import pandas as pd

if __name__ == '__main__':
    metadata = pd.read_table(sys.argv[1]).set_index("sample_id")
    grouping_column = sys.argv[2]


    metadata[grouping_column] = metadata[grouping_column].astype(str)

    all_terms = metadata[grouping_column].unique()
    pd.DataFrame(
        {
            "value": all_terms,
            "path_safe_id": [
                str(i) for i in range(len(all_terms))
            ],
            "column_name": [grouping_column] * len(all_terms)
        }
    ).to_csv(
        sys.argv[3],
        sep='\t',
        index=False
    )