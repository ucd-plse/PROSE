# %%
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from glob import glob
import xarray as xr
from copy import deepcopy
import subprocess
import pickle

# %%
VAR_NAMES = []
def plot_performance_vs_correctness_L2(df, error_type, error_threshold, optimal=[]):
    
    df_plot = df.dropna(subset=["Cost", f"{error_type}, Relative Error (L2 Norm)"])

    baseline_cost = df[(df['Configuration Number'] == 0)]['Cost'].iloc[0]
    df_plot = df_plot.assign(Improvement = baseline_cost/df_plot['Cost'])
    
    fig = px.scatter(
        df_plot,
        x = "Improvement",
        y = f"{error_type}, Relative Error (L2 Norm)",
        color = '32-bit %',
        hover_data = ['Configuration Number'],
        log_y=True,
    )

    fig.update_traces(
        marker={
            'size' : 14,
            'opacity' : 0.5,
        }
    )


    fig.update_traces(
        marker={
            'size' : 14,
            'line_width' : 1,
            'line_color' : "black",
        }
    )
    
    if not df_plot[df_plot["Configuration Number"] == 1].empty:
        fig.add_trace(
            go.Scatter(
                x = df_plot[df_plot["Configuration Number"] == 1]["Improvement"],
                y = df_plot[df_plot["Configuration Number"] == 1][f"{error_type}, Relative Error (L2 Norm)"],
                mode = "markers+text",
                marker_symbol = 'circle-open',
                marker_size = 14,
                marker_color = "black",
                marker_line_width = 3,
                marker_line_color = "black",
                name = "Uniform 32-bit",
                showlegend=False,
                text=["Uniform 32-bit"],
                textposition=["middle left"],
            )
        )

    fig.update_layout(
        title=dict(
            text = "MPAS-A",# (Entire Model)",
        ),
        yaxis_tickformat = ".1e",
        xaxis_tickformat = ".1f",
        xaxis_ticksuffix = "x",
        xaxis_title = "Speedup",
        yaxis_title = "Relative Error",
        font_family = "Times New Roman",

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

    return fig


def get_MPAS_data(search_log_path):

    global LONG_NAME_MAP
    global VAR_NAMES

    with open(search_log_path, "r") as f:
        search_log_lines = f.readlines()

    df_entire = []
    df_subset = []
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
        config = []
        VAR_NAMES = []
        with open(row['Configuration Path'], "r" ) as f:
            llines = f.readlines()
            for lline in llines:
                if ",4" in lline:
                    float_count += 1
                    config.append(0)
                elif ",8" in lline:
                    double_count += 1
                    config.append(1)

                VAR_NAMES.append(lline.split(",")[0])

        row["config"] = config
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

        try:
            errors_df = pd.read_pickle(os.path.join(config_dir_path, "errors.pckl"))
            for column_name in errors_df.columns:
                if column_name != "time":
                    # row[f"{metric}, Relative Error (Average)"] = np.mean(errors_df[column_name])
                    # row[f"{metric}, Relative Error (Variance)"] = np.var(errors_df[column_name])
                    # row[f"{metric}, Relative Error (Median)"] = np.median(errors_df[column_name])
                    # row[f"{metric}, Relative Error (Max)"] = np.max(errors_df[column_name])
                    # row[f"{metric}, Relative Error (Min)"] = np.min(errors_df[column_name])
                    # row[f"{metric}, Relative Error (75th percentile)"] = np.percentile(errors_df[column_name], 75)
                    # row[f"{metric}, Relative Error (25th percentile)"] = np.percentile(errors_df[column_name], 25)
                    row[f"{column_name}, Relative Error (L2 Norm)"] = np.linalg.norm(errors_df[column_name], ord=2)

        except FileNotFoundError:
            for column_name in errors_df.columns:
                if column_name != "time":
                    # row[f"{metric}, Relative Error (Average)"] = np.nan
                    # row[f"{metric}, Relative Error (Variance)"] = np.nan
                    # row[f"{metric}, Relative Error (Median)"] = np.nan
                    # row[f"{metric}, Relative Error (Max)"] = np.nan
                    # row[f"{metric}, Relative Error (Min)"] = np.nan
                    # row[f"{metric}, Relative Error (75th percentile)"] = np.nan
                    # row[f"{metric}, Relative Error (25th percentile)"] = np.nan
                    row[f"{column_name}, Relative Error (L2 Norm)"] = np.nan

        try:
            row['Cost'] = float(line.split()[4])
        except ValueError:
            try:
                row['Cost'] = float(line.split()[5])
            except ValueError:
                row['Cost'] = np.nan

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
                        if lline.strip().lower().startswith(key.lower() + "::") or ("fluxes" in row['Procedure Name'] and "flux" in lline.lower()):
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
# LONG_NAME_MAP = xr.open_dataset("./prose_logs/0000/history.2014-09-10_00.00.00.nc")
df_entire, df_subset = get_MPAS_data("./prose_logs/search_log.txt")

# %%
fig = plot_performance_vs_correctness_L2(df_entire, error_type="ke", error_threshold=132)
fig.write_html("figure_5_MPAS.html")

# %%
from plotly.subplots import make_subplots

df_subset["Average CPU Time Per Call"] = df_subset["CPU Time"] / df_subset["Execution Count"]
for procedure_name in df_subset[df_subset["Procedure Name"].str.contains("::")]["Procedure Name"].unique():
    for config_hash in df_subset[df_subset["Procedure Name"] == procedure_name]["config_hash"].unique():
        df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"] = df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"].mean()

procedure_percentages = {
    '::atm_time_integration::atm_recover_large_step_variables_work': "23.9",#23.85573770006054,
    '::atm_time_integration::atm_compute_dyn_tend_work': "21.5",#21.513575486114263,
    '::atm_time_integration::atm_advance_acoustic_step_work': "9.9",#9.919105842657316,
    '::atm_time_integration::fluxes': "2.6",#2.598569413169517,
}

subplot_titles = [f'{x[x.rfind(":") + 1:]} ({procedure_percentages[x]}%)' if "fluxes" not in x else f'flux3 and flux4 ({procedure_percentages[x]}%)' for x in procedure_percentages.keys()]

fig = make_subplots(len(procedure_percentages),1, subplot_titles=tuple(subplot_titles), shared_xaxes=True, vertical_spacing=0.1)
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
        text = "MPAS-A",
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

fig.update_xaxes(type="log",ticksuffix="x")
fig.update_xaxes(title="Speedup", row=len(procedure_percentages), col=1)
fig.update_yaxes(visible=False)
fig.update_annotations(yshift=-5, font_size=18)
fig.write_html("figure_6_MPAS.html")

# %%



