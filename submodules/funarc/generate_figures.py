# %%
import pandas as pd
import numpy as np
import plotly.express as px
import os
from glob import glob

# %%
def plot_performance_vs_correctness(df):
    
    df_plot = df.dropna()
    
    baseline_cost = df[(df['Configuration Number'] == 0)]['Cost'].iloc[0]
    baseline_output = df[(df['Configuration Number'] == 0)]['Output'].iloc[0]
    df_plot = df_plot.assign(Improvement = baseline_cost/df_plot['Cost'])
    df_plot = df_plot.assign(Error = ((baseline_output - df_plot['Output'])/baseline_output).abs())

    fig = px.scatter(
        df_plot,
        x="Improvement",
        y="Error",
        hover_data = ['Configuration Number'],
        color="32-bit %"
    )
    fig.update_traces(
        marker={
            'size' : 14,
            'opacity' : 0.5,
            'line_width' : 1,
            'line_color' : "black",
        },
        showlegend = False,
    )

    fig.add_vline(x=1.0, line_width=2, line_dash="dash", line_color="grey")
    fig.add_hline(y=3e-4, line_width=2, line_dash="dash", line_color="grey")

    fig.update_layout(
        yaxis_tickformat = ".0e",
        xaxis_tickformat = ".1f",
        xaxis_ticksuffix = "x",
        xaxis_title = "Speedup",
        yaxis_title = "Relative Error",
        font_family = "Times New Roman",
        legend = dict(
            title_text = ""
        ),
    coloraxis={
        "cmin" : 0,
        "cmax" : 100,
        "colorbar" : {
            'title' : "% 32-bit",
        },
    },
    )


    return fig


def plot_label_frequency(df):

    fig = px.histogram(
        df,
        x = "Label",
    )
    fig.update_layout(
        font_size = 18
    )
    
    return fig


def get_funarc_data(search_log_path):
    
    with open(search_log_path, "r") as f:
        search_log_lines = f.readlines()

    df = []
    for line in search_log_lines:
        row = {}

        try:
            row['Configuration Number'] = int(line.split(":")[0])
        except ValueError:
            continue
        config_dir_path = os.path.join(os.path.dirname(search_log_path), f"{row['Configuration Number']:0>4}")
        config_path = glob(f"{config_dir_path}/config*")
        assert ( len(config_path) == 1)
        row['Configuration Path'] = config_path[0]

        float_count = 0
        double_count = 0
        with open(row['Configuration Path'], "r" ) as f:
            llines = f.readlines()
            for lline in llines:
                if ",4" in lline:
                    float_count += 1
                elif ",8" in lline:
                    double_count += 1

        row['32-bit %'] = 100*(float_count / (double_count + float_count))

        if "[PASSED]" in line:
            row['Label'] = "Passed"
            row['Cost'] = float(line.split()[4])
        elif "error threshold was exceeded" in line:
            row['Label'] = "Failed"
            row['Cost'] = float(line.split()[4])
        elif "(timeout)" in line:
            row['Label'] = "Timeout"
            row['Cost'] = float(line.split()[4])
        elif "(runtime failure)" in line:
            row['Label'] = "Runtime Error"
            row['Cost'] = np.nan
        elif "(compilation error)" in line:
            row['Label'] = 'Compilation Error'
            row['Cost'] = np.nan
        elif "(plugin error)" in line:
            row['Label'] = 'Prose Plugin Error'
            row['Cost'] = np.nan
        else:
            continue

        try:
            with open(os.path.join(config_dir_path,"outlog.txt"), "r") as f:
                
                line = f.readlines()[0]
                if line.startswith(" out:"):
                    row["Output"] = float(line.strip().split()[-1])
                else:
                    row["Output"] = np.nan

        except (FileNotFoundError, IndexError):
            row["Output"] = np.nan

        df.append(row)

    return pd.DataFrame(df)

# %%
df = get_funarc_data("./prose_logs/search_log.txt")

# %%
fig = plot_performance_vs_correctness(df)
fig.write_html("funarc_search_results.html")

# %%



