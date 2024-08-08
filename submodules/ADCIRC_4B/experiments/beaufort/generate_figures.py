# %%
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from glob import glob
import pickle
from copy import deepcopy
import subprocess

# %%
def plot_performance_vs_correctness(df, error_type, error_threshold):
    
    df_plot = df.dropna()

    baseline_total_cost = df_plot[(df['Configuration Number'] == 0)]['Total Cost'].iloc[0]
    baseline_subset_cost = df_plot[(df_plot['Configuration Number'] == 0)]['Subset Cost'].iloc[0]
    df_plot = df_plot.assign(Improvement = baseline_subset_cost/df_plot['Subset Cost'])
    df_plot = df_plot.assign(Total_Improvement = baseline_total_cost/df_plot['Total Cost'])
    df_plot = df_plot.assign(error_x_minus = df_plot['Improvement'] - df_plot['Total_Improvement'])
    df_plot = df_plot.assign(error_x = df_plot['Total_Improvement'] - df_plot['Improvement'])

    num = df_plot._get_numeric_data()
    num[num < 0] = 0

    fig = px.scatter(
        df_plot,
        x="Improvement",
        y=error_type,
        color = '32-bit %',
        hover_data = ['Configuration Number'],
        log_y = True, 
    )
    fig.update_traces(
        marker={
            'size' : 14,
            'opacity' : 0.5,
            'line_width' : 1,
            'line_color' : "black",
        }
    )
    if not df_plot[df_plot["Configuration Number"] == 1].empty:
        fig.add_trace(
            go.Scatter(
                x = df_plot[df_plot["Configuration Number"] == 1]["Improvement"],
                y = df_plot[df_plot["Configuration Number"] == 1][error_type],
                mode = "markers+text",
                marker_symbol = 'circle-open',
                marker_size = 14,
                marker_color = "black",
                marker_line_width = 3,
                marker_line_color = "black",
                name = "Uniform 32-bit",
                showlegend=False,
                text=["Uniform 32-bit"],
                textposition=["top right"],
            )
        )

    fig.update_layout(
        title=dict(
            text = "ADCIRC",
        ),
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
            'title' : "% 32-bit (Hotspot)",
        },
    },
    )
    fig.add_vline(x=1.0, line_width=2, line_dash="dash", line_color="grey")
    fig.add_hline(y=error_threshold, line_width=2, line_dash="dash", line_color="grey")

    return fig, df_plot

def get_ADCIRC_data(search_log_path):
    
    with open(search_log_path, "r") as f:
        search_log_lines = f.readlines()

    df_entire = []
    df_subset = []
    for line in search_log_lines:
        print(line)
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
        elif "error threshold was exceeded" in line:
            row['Label'] = "Exceeded Error Threshold"
        elif "(timeout)" in line:
            row['Label'] = "Timeout"
        elif "(runtime failure)" in line:
            row['Label'] = "Runtime Error"
        elif "(compilation error)" in line:
            row['Label'] = 'Compilation Error'
        elif "(plugin error)" in line:
            row['Label'] = 'Prose Plugin Error'
        else:
            continue

        if row["Configuration Number"] == 0:
            row["Water Elevation Error (L2 Norm)"] = 0
            row["Water Elevation Error (Max Norm)"] = 0
            row["Water Velocity Error (L2 Norm)"] = 0
            row["Water Velocity Error (Max Norm)"] = 0
        else:
            try:
                with open(os.path.join(config_dir_path,"error_metrics.txt"), "r") as f:
                    temp = f.readlines()[0]
                    temp = temp.split()
                    assert( len(temp) == 7 )
                    row["Water Elevation Error (L2 Norm)"] = float(temp[3])
                    row["Water Elevation Error (Max Norm)"] = float(temp[4])
                    row["Water Velocity Error (L2 Norm)"] = float(temp[5])
                    row["Water Velocity Error (Max Norm)"] = float(temp[6])

            except FileNotFoundError:
                row["Water Elevation Error (L2 Norm)"] = np.nan
                row["Water Elevation Error (Max Norm)"] = np.nan
                row["Water Velocity Error (L2 Norm)"] = np.nan
                row["Water Velocity Error (Max Norm)"] = np.nan

        tokens = [token.strip() for token in line.split()]
        try:
            row['Subset Cost'] = float(tokens[tokens.index("=") + 1])
            row['Total Cost'] = float(tokens[tokens.index("=", tokens.index("=") + 1) + 1][:-1])
        except ValueError:
            row['Subset Cost'] = np.nan
            row['Total Cost'] = np.nan

        df_entire.append(deepcopy(row))

        try:
            with open(os.path.join(config_dir_path, "gptl_subset_info.pckl"), "rb") as f:
                gptl_subset_info = pickle.load(f)
                
            subprocess.run(
                f"tar -xvf gptl_timing.tar.gz",
                shell = True,
                cwd = config_dir_path,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            execution_counts = {}
            with open(os.path.join(config_dir_path, "timing.000000"), "r") as f:
                for line in f.readlines():
                    if line[2:].lstrip().startswith("::"):
                        procedure_name = line[2:].split()[0].strip().lower()
                        execution_count = float(line[2:].split()[1].strip())
                        execution_counts[procedure_name] = execution_count

            subprocess.run(
                f"rm timing.*",
                shell = True,
                cwd = config_dir_path,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            for key, value in gptl_subset_info.items():
                row['Procedure Name'] = key
                row['CPU Time'] = value
                try:
                    row["Execution Count"] = execution_counts[row['Procedure Name'].lower()]
                except KeyError:
                    row["Execution Count"] = np.nan
                float_count = 0
                double_count = 0
                config_hash = ""
                with open(row['Configuration Path'], "r" ) as f:
                    llines = f.readlines()
                    for lline in llines:
                        if lline.strip().lower().startswith(key.lower() + "::"):
                            if ",4" in lline:
                                float_count += 1
                                config_hash = config_hash + "4"
                            elif ",8" in lline:
                                double_count += 1
                                config_hash = config_hash + "8"
                
                row["config_hash"] = hash(config_hash)
                try:
                    row['32-bit %'] = 100*(float_count / (double_count + float_count))
                except Exception as e:
                    row['32-bit %'] = np.nan

                df_subset.append(deepcopy(row))

        except FileNotFoundError:
            continue

    return pd.DataFrame(df_entire), pd.DataFrame(df_subset)

# %%
df_entire, df_subset = get_ADCIRC_data("./prose_logs/search_log.txt")

# %%
fig, _ = plot_performance_vs_correctness(df_entire, error_type="Water Elevation Error (L2 Norm)", error_threshold=1e-1)
fig.write_html("figure_5_ADCIRC.html")

# %%
from plotly.subplots import make_subplots

df_subset["Average CPU Time Per Call"] = df_subset["CPU Time"] / df_subset["Execution Count"]
for procedure_name in df_subset[df_subset["Procedure Name"].str.contains("::")]["Procedure Name"].unique():
    for config_hash in df_subset[df_subset["Procedure Name"] == procedure_name]["config_hash"].unique():
        df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"] = df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"].mean()

procedure_percentages = {
 '::itpackv::peror': "36.3", #36.29100455739751,
 '::itpackv::pjac': "32.9", #31.977778859505513,
 '::itpackv::jcg': "<0.01",#2.286905914467208e-05,
}

fig = make_subplots(len(procedure_percentages),1, subplot_titles=tuple([f'{x[x.rfind(":") + 1:]} ({procedure_percentages[x]}%)' for x in procedure_percentages.keys()]), shared_xaxes=True, vertical_spacing=0.15)
for i, procedure_name in enumerate(procedure_percentages.keys()):
    df_plot = df_subset[df_subset["Procedure Name"] == procedure_name]
    baseline_cost = df_plot[df_plot["Configuration Number"] == 0]["Average CPU Time Per Call"].values[0]
    df_plot = df_plot.assign(Improvement = np.round(baseline_cost/df_plot['Average CPU Time Per Call'], decimals=2))
    df_plot = df_plot.drop_duplicates(subset=["config_hash","Improvement"])
    fig.add_trace(
        go.Scatter(
            x = df_plot["Improvement"],
            y = np.random.rand(len(df_subset)),
            mode = 'markers',
            customdata=df_plot["Configuration Number"],
            hovertemplate="%{customdata}",
            marker = dict(
                size = 10,
                color=df_plot["32-bit %"],
                line_width = 1,
                line_color = "black",
                opacity = 0.6,
                coloraxis="coloraxis1",
                symbol = "diamond",
            ),
            showlegend=False
        ),
        i + 1,
        1
    )
    if i == 0:
        fig.update_traces(
            marker_colorbar_title = "% 32-bit",
            marker_colorscale = "Plasma",
        )

fig.update_layout(
    showlegend = False,
    title=dict(
        text = "ADCIRC",
    ),
    font_family = "Times New Roman",
    coloraxis={
        "cmin" : 0,
        "cmax" : 100,
        "colorbar" : {
            'title' : "% 32-bit (Procedure)",
        },
    },
)

fig.update_xaxes(type="log", ticksuffix="x")
fig.update_xaxes(title="Speedup", row=len(procedure_percentages), col=1)
fig.update_yaxes(visible=False)
fig.update_annotations(yshift=-5, font_size=18)
fig.write_html("figure_6_ADCIRC.html")