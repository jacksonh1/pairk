import pandas as pd

def filter_critical_errors_out(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["critical_error"].isna()].copy()


