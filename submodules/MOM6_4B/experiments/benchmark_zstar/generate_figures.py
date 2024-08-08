# %%
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from glob import glob
import pickle
from copy import deepcopy

# %%
ERROR_METRICS = [
    "Energy/Mass",
    "Maximum CFL",
    "Mean Sea Level",
    "Total Mass",
    "Mean Salin",
    "Frac Mass Err",
    "Salin Err",
    "Temp Err",
]

def plot_performance_vs_correctness(df, error_type, norm):
    
    df_plot = df.dropna(subset=["Cost", f"{error_type} ({norm})"])
    baseline_cost = df_plot[(df_plot['Configuration Number'] == 0)]['Cost'].iloc[0]
    df_plot = df_plot.assign(Improvement = baseline_cost/df_plot['Cost'])
    
    fig = px.scatter(
        df_plot,
        x="Improvement",
        y=f"{error_type} ({norm})",
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
                y = df_plot[df_plot["Configuration Number"] == 1][f"{error_type} ({norm})"],
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
            text = "MOM6",
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
    fig.add_hline(y=0.25, line_width=2, line_dash="dash", line_color="grey")

    return fig, df_plot


def plot_label_frequency(df):

    fig = px.histogram(
        df,
        x = "Label",
    )
    fig.update_layout(
        font_size = 18
    )

    return fig


def get_MOM6_data(search_log_path):

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

        try:
            with open(os.path.join(config_dir_path, "errors_L2.npy"), "rb") as f:
                errors = np.load(f)
            for i, metric in enumerate(ERROR_METRICS):
                row[f"{metric} (L2 norm)"] = errors[i]

            with open(os.path.join(config_dir_path, "errors_inf.npy"), "rb") as f:
                errors = np.load(f)
            for i, metric in enumerate(ERROR_METRICS):
                row[f"{metric} (inf norm)"] = errors[i]

        except FileNotFoundError:
            for metric in ERROR_METRICS:
                row[f"{metric} (L2 norm)"] = np.nan
                row[f"{metric} (inf norm)"] = np.nan
 
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

            import subprocess
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
df_entire, df_subset = get_MOM6_data("./prose_logs/search_log.txt")

# %%
fig, df_plot = plot_performance_vs_correctness(df_entire, error_type=ERROR_METRICS[1], norm="L2 norm")
fig.write_html("figure_5_MOM6.html")

# %%
from plotly.subplots import make_subplots

df_subset["Average CPU Time Per Call"] = df_subset["CPU Time"] / df_subset["Execution Count"]
for procedure_name in df_subset[df_subset["Procedure Name"].str.contains("::")]["Procedure Name"].unique():
    if "flux_adjust" in procedure_name:
        continue # don't take the average of all 64-bit variants of flux-adjust routines which are noisy
    for config_hash in df_subset[df_subset["Procedure Name"] == procedure_name]["config_hash"].unique():
        df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"] = df_subset.loc[(df_subset["Procedure Name"] == procedure_name) & (df_subset["config_hash"] == config_hash), "Average CPU Time Per Call"].mean()

procedure_percentages = {
    '::mom_continuity_ppm::zonal_flux_layer': "22.4",#22.43212079470703,
    '::mom_continuity_ppm::zonal_mass_flux': "9.0",#8.9596406665791,
    '::mom_continuity_ppm::zonal_flux_adjust': "7.6",#7.572937598608928,
    '::mom_continuity_ppm::set_zonal_bt_cont': "6.3",#6.251864356826961,
    '::mom_continuity_ppm::zonal_face_thickness': "2.3",#2.3120783396943727,
}

fig = make_subplots(len(procedure_percentages),1, subplot_titles=tuple([f'{x[x.rfind(":") + 1:]} ({procedure_percentages[x]}%)' for x in procedure_percentages.keys()]), shared_xaxes=True, vertical_spacing=0.09)
for i, procedure_name in enumerate(procedure_percentages.keys()):
    df_plot = df_subset[df_subset["Procedure Name"] == procedure_name]
    baseline_cost = df_plot[df_plot["Configuration Number"] == 0]["Average CPU Time Per Call"].values[0]
    df_plot = df_plot.assign(Improvement = np.round(baseline_cost/df_plot['Average CPU Time Per Call'], decimals=2))
    if "flux_adjust" in procedure_name: # reduce crowding of flux adjust variants with > 0.1x improvement
        df_plot.loc[df_plot["Improvement"] > 0.1, "Improvement"] = np.round(df_plot[df_plot["Improvement"] > 0.1]["Improvement"], decimals=1)
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
        text = "MOM6",
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
fig.write_html("figure_6_MOM6.html")
