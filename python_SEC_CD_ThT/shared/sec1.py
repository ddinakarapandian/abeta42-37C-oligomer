import os
import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scipy.signal import find_peaks
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

class SECData:
    """
    A class for handling and analyzing Size Exclusion Chromatography (SEC) data.
    """

    def __init__(self, data):
        """
        Initialize the SECData object.

        Args:
            data (list, pd.DataFrame): Input data, can be a list of DataFrames, 
                                     a list of filenames, or a single DataFrame.
        """
        self.data = self._initialize_data(data)
        
    def _initialize_data(self, data):
        """Helper method to initialize data with proper error handling."""
        if isinstance(data, pd.DataFrame):
            return [data]
        elif isinstance(data, list):
            if all(isinstance(item, pd.DataFrame) for item in data):
                return data
            elif all(isinstance(item, str) for item in data):
                return self._load_data_from_files(data)
            else:
                raise ValueError("List must contain only DataFrames or only filenames")
        elif isinstance(data, str):
            return self._load_data_from_files([data])
        else:
            raise ValueError("Input must be a DataFrame, a list of DataFrames, a list of filenames, or a single filename")
    
    def _load_data_from_files(self, filenames):
        """Helper method to load data from files with flexible parsing."""
        loaded_data = []
        for filename in filenames:
            df = None
            # Try different loading strategies
            for encoding in ['utf-16', 'utf-8']:
                try:
                    # First try tab-separated
                    df = pd.read_csv(filename, encoding=encoding, sep='\t', skiprows=3)
                    if len(df.columns) >= 2:  # If successful, break
                        break
                        
                    # If only one column, try comma-separated
                    df = pd.read_csv(filename, encoding=encoding, skiprows=3)
                    if len(df.columns) >= 2:  # If successful, break
                        break
                        
                    # If still one column, check if values are comma-separated within the column
                    if len(df.columns) == 1:
                        # Try splitting the first column on commas
                        first_row = df.iloc[0, 0]
                        if isinstance(first_row, str) and ',' in first_row:
                            # Split the single column into multiple columns
                            new_df = df[df.columns[0]].str.split(',', expand=True)
                            # Convert to numeric, coercing errors to NaN
                            new_df = new_df.apply(pd.to_numeric, errors='coerce')
                            df = new_df
                            break
                except Exception as e:
                    continue
                    
            if df is None:
                raise ValueError(f"Failed to load file {filename} with any parsing strategy")
                
            # Ensure we have at least two numeric columns
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 2:
                raise ValueError(f"File {filename} does not contain at least two numeric columns")
                
            # Take the first two numeric columns
            df = df[numeric_cols[:2]]
            loaded_data.append(df)
            
        return loaded_data

    def rename_columns(self):
        """Rename columns of all DataFrames in self.data."""
        for i, df in enumerate(self.data):
            try:
                # Ensure we have numeric data in at least two columns
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) < 2:
                    raise ValueError(f"DataFrame {i} does not have enough numeric columns")
                
                # Take first two numeric columns and rename them
                new_df = df[numeric_cols[:2]].copy()
                new_df.columns = ['Volume_ml', 'Absorbance_mAU']
                self.data[i] = new_df
            except Exception as e:
                raise ValueError(f"Error processing DataFrame {i}: {str(e)}")

    def filter_by_volume(self, min_volume=0, max_volume=25):
        """Filter data by volume range."""
        filtered_data = [df[(df['Volume_ml'] >= min_volume) & (df['Volume_ml'] <= max_volume)] 
                        for df in self.data]
        return SECData(filtered_data)

    def normalize_data(self):
        """Normalize absorbance data for all DataFrames."""
        for df in self.data:
            min_val = df['Absorbance_mAU'].min()
            max_val = df['Absorbance_mAU'].max()
            df['Normalized_Absorbance_mAU'] = (df['Absorbance_mAU'] - min_val) / (max_val - min_val)
        return self

    def adjust_baseline(self, method='simple', baseline_regions=None):
        """
        Adjust baseline using either simple minimum subtraction or piecewise correction.
        Handles both positive slope and negative offsets in baseline.
        
        Args:
            method (str): Either 'simple', 'linear', or 'piecewise'
            baseline_regions (list): List of tuples (start_vol, end_vol) for piecewise correction
            
        Returns:
            SECData: Self with adjusted data
        """
        if method == 'simple':
            # Simple offset correction
            for df in self.data:
                min_value = df['Absorbance_mAU'].min()
                df['Adjusted_Absorbance_mAU'] = df['Absorbance_mAU'] - min_value
                
        elif method == 'linear':
            # Linear trend correction using start and end points
            for df in self.data:
                # Use first and last 5% of points to determine baseline trend
                n_points = len(df)
                n_edge = max(int(n_points * 0.1), 50)  # At least 5 points
                
                start_points = df.iloc[:n_edge]
                end_points = df.iloc[-n_edge:]
                
                # Fit line through start and end regions
                x = np.concatenate([start_points['Volume_ml'], end_points['Volume_ml']])
                y = np.concatenate([start_points['Absorbance_mAU'], end_points['Absorbance_mAU']])
                coeffs = np.polyfit(x, y, 1)
                baseline = np.polyval(coeffs, df['Volume_ml'])
                
                # Subtract baseline
                df['Adjusted_Absorbance_mAU'] = df['Absorbance_mAU'] - baseline
                
        elif method == 'piecewise':
            if baseline_regions is None:
                baseline_regions = [
                    (0, 6),     # Pre-peak region
                    (12, 13),    # Between first two peaks
                    (19, 20),   # Between second and third peaks
                    (22, 25)    # Post-peak region
                ]
            
            adjusted_data = []
            for df in self.data:
                adjusted_df = df.copy()
                
                # First find global baseline trend
                all_baseline_points = pd.DataFrame()
                for start_vol, end_vol in baseline_regions:
                    mask = (adjusted_df['Volume_ml'] >= start_vol) & (adjusted_df['Volume_ml'] <= end_vol)
                    all_baseline_points = pd.concat([all_baseline_points, adjusted_df[mask]])
                
                # Fit global trend if we have enough baseline points
                if len(all_baseline_points) > 2:
                    global_coeffs = np.polyfit(all_baseline_points['Volume_ml'], 
                                            all_baseline_points['Absorbance_mAU'], 1)
                    global_baseline = np.polyval(global_coeffs, adjusted_df['Volume_ml'])
                    adjusted_df['Absorbance_mAU'] -= global_baseline
                
                # Then do local corrections
                for start_vol, end_vol in baseline_regions:
                    mask = (adjusted_df['Volume_ml'] >= start_vol) & (adjusted_df['Volume_ml'] <= end_vol)
                    region_data = adjusted_df[mask]
                    
                    if not region_data.empty:
                        # Local baseline correction
                        x = region_data['Volume_ml'].values
                        y = region_data['Absorbance_mAU'].values
                        local_coeffs = np.polyfit(x, y, 1)
                        local_baseline = np.polyval(local_coeffs, adjusted_df['Volume_ml'])
                        
                        # Apply local correction after this region
                        mask_after = adjusted_df['Volume_ml'] >= start_vol
                        adjusted_df.loc[mask_after, 'Absorbance_mAU'] -= local_baseline[mask_after]
                
                adjusted_data.append(adjusted_df)
            
            self.data = adjusted_data
            
        else:
            raise ValueError("Method must be 'simple', 'linear', or 'piecewise'")
            
        return self
        
    def export_data(self, output_filename, include_placeholder=True):
        """
        Export the current data to a CSV file.
        
        Args:
            output_filename (str): Name of the output file
            include_placeholder (bool): Whether to add placeholder rows at the beginning
            
        Returns:
            pd.DataFrame: The exported dataframe
        """
        # Combine all dataframes
        combined_df = pd.concat(self.data, ignore_index=True)
        
        if include_placeholder:
            # Add placeholder rows at the beginning
            placeholder_rows = pd.DataFrame({col: ['Placeholder'] * 2 for col in combined_df.columns})
            combined_df = pd.concat([placeholder_rows, combined_df], ignore_index=True)
        
        # Ensure Volume_ml and Absorbance_mAU are the first two columns
        first_columns = ['Volume_ml', 'Absorbance_mAU']
        other_columns = [col for col in combined_df.columns if col not in first_columns]
        combined_df = combined_df[first_columns + other_columns]
        
        # Export to CSV
        combined_df.to_csv(output_filename, index=False)
        print(f"Data exported to {output_filename}")
        
        return combined_df

    def process_data_to_arrays(self, min_volume=0, max_volume=25, normalize=True):
        """Process data by filtering, normalizing, and converting to arrays."""
        processed_arrays = []
        for df in self.data:
            # Ensure proper column names
            if len(df.columns) < 2:
                raise ValueError("DataFrame has insufficient columns")
                
            df = df.copy()
            df.columns = ['Volume_ml', 'Absorbance_mAU'] + [f'Unnamed_{j}' for j in range(len(df.columns) - 2)]
            
            # Filter and process
            df = df[(df['Volume_ml'] >= min_volume) & (df['Volume_ml'] <= max_volume)]
            self._process_single_dataframe(df, normalize)
            
            # Select appropriate column for output
            y_col = 'Normalized_Absorbance_mAU' if normalize else 'Adjusted_Absorbance_mAU'
            processed_arrays.append(df[['Volume_ml', y_col]].to_numpy())
            
        return processed_arrays

    def _process_single_dataframe(self, df, normalize):
        """Helper method to process a single DataFrame."""
        # Baseline adjustment
        min_value = df['Absorbance_mAU'].min()
        df['Adjusted_Absorbance_mAU'] = df['Absorbance_mAU'] - min_value if min_value < 0 else df['Absorbance_mAU']
        
        # Normalization if requested
        if normalize:
            min_val = df['Adjusted_Absorbance_mAU'].min()
            max_val = df['Adjusted_Absorbance_mAU'].max()
            df['Normalized_Absorbance_mAU'] = (df['Adjusted_Absorbance_mAU'] - min_val) / (max_val - min_val)

    def plot(self, data=None, normalized=False, title='', xaxis_title='', yaxis_title='', 
            legend_title='', x_range=None, y_range=None, trace_names=None, legend=None):
        """
        Create an interactive Plotly plot of the SEC data.
        
        Args:
            data (SECData, optional): Data to plot. If None, uses self.data
            normalized (bool): Whether to use normalized data
            title (str): Plot title
            xaxis_title (str): X-axis title
            yaxis_title (str): Y-axis title
            legend_title (str): Legend title
            x_range (list): X-axis range [min, max]
            y_range (list): Y-axis range [min, max]
            trace_names (list): Names for each trace
            legend (dict): Legend parameters
            
        Returns:
            plotly.graph_objs.Figure: The created figure
        """
        fig = go.Figure()
        
        # Determine which data to plot
        if data is None:
            data_to_plot = self.data
        elif isinstance(data, SECData):
            data_to_plot = data.data
        else:
            data_to_plot = [data] if isinstance(data, pd.DataFrame) else data
            
        # Choose column for y-axis
        y_column = 'Normalized_Absorbance_mAU' if normalized else 'Absorbance_mAU'
        
        # Generate trace names if not provided
        trace_names = trace_names or [f'Data {i+1}' for i in range(len(data_to_plot))]
        
        # Add each trace to the figure
        for name, df in zip(trace_names, data_to_plot):
            try:
                fig.add_trace(go.Scatter(
                    x=df['Volume_ml'],
                    y=df[y_column],
                    mode='lines',
                    name=name
                ))
            except KeyError as e:
                print(f"Warning: Could not find column {str(e)}. Available columns: {df.columns.tolist()}")
                continue
                
        # Update layout
        layout_updates = {
            'title': title,
            'xaxis_title': xaxis_title or 'Volume (ml)',
            'yaxis_title': yaxis_title or ('Normalized Absorbance' if normalized else 'Absorbance (mAU)'),
            'legend_title': legend_title,
            'showlegend': legend is not None
        }
        
        if x_range:
            layout_updates['xaxis'] = {'range': x_range}
        if y_range:
            layout_updates['yaxis'] = {'range': y_range}
            
        # Set default legend position if not provided
        if legend is None:
            legend = dict(yanchor="top", y=0.99, xanchor="right", x=0.99)
        layout_updates['legend'] = legend
        
        fig.update_layout(**layout_updates)
        
        # Display the figure
        fig.show()
        
    def _update_plot_layout(self, fig, title, xaxis_title, yaxis_title, legend_title, x_range, y_range, legend):
        """Helper method to update plot layout."""
        if x_range:
            fig.update_xaxes(range=x_range)
        if y_range:
            fig.update_yaxes(range=y_range)
            
        legend = legend or dict(yanchor="top", y=0.99, xanchor="right", x=0.99)
        
        fig.update_layout(
            title=title,
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title,
            legend_title=legend_title,
            legend=legend
        )
    def list_plot(self, data=None, normalized=False, title='', xlabel='', ylabel='', 
                  legend_labels=None, colors=None, linestyles=None, linewidth=2.5, linewidths=None,
                  x_range=None, y_range=None, figsize=(12, 10), dpi=300, output_file=None,
                  x_major_ticks=None, y_major_ticks=None, x_minor_ticks=None, y_minor_ticks=None,
                  major_tick_length=10, minor_tick_length=6, major_tick_width=2, minor_tick_width=1,
                  show_box=True, legend_box=False,
                  legend_position=None, legend_fontsize=24,
                  title_fontsize=24, title_weight='bold', title_font='Arial',
                  xlabel_fontsize=24, xlabel_weight='normal', xlabel_font='Arial',
                  ylabel_fontsize=24, ylabel_weight='normal', ylabel_font='Arial'):
        """Create a publication-quality Matplotlib plot."""
        plt.rcParams['font.family'] = 'Arial'
        
        # Process data
        if data is None:
            data = self.process_data_to_arrays(normalize=normalized)
        elif isinstance(data, SECData):
            data = data.process_data_to_arrays(normalize=normalized)
        data = [data] if not isinstance(data, list) else data
        
        # Create figure and plot
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        self._setup_plot_style(ax, data, colors, linestyles, linewidth, linewidths, legend_labels)
        self._configure_axes(ax, normalized, x_major_ticks, y_major_ticks, x_minor_ticks, y_minor_ticks,
                           major_tick_length, minor_tick_length, major_tick_width, minor_tick_width, show_box)
        self._set_plot_labels(ax, title, xlabel, ylabel, normalized, title_fontsize, title_weight, title_font, xlabel_fontsize, xlabel_weight, xlabel_font, ylabel_fontsize, ylabel_weight, ylabel_font)
        
        # Set ranges and finalize
        if x_range:
            ax.set_xlim(x_range)
        if y_range:
            ax.set_ylim(y_range)
            
        # Configure legend
        if legend_position:
            ax.legend(fontsize=legend_fontsize, frameon=legend_box, bbox_to_anchor=legend_position)
        else:
            ax.legend(fontsize=legend_fontsize, frameon=legend_box)
            
        # Make legend labels bold
        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_weight('bold')
            
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
            print(f"Figure saved as {output_file}")
        else:
            plt.show()


    def _setup_plot_style(self, ax, data, colors, linestyles, linewidth, linewidths, legend_labels):
        """Configure plot style parameters."""
        num_datasets = len(data)
        
        # Ensure all style parameters have the correct length
        if colors is not None:
            colors = colors * (num_datasets // len(colors) + 1)
            colors = colors[:num_datasets]
        else:
            colors = plt.cm.viridis(np.linspace(0, 1, num_datasets))
            
        if linestyles is not None:
            linestyles = linestyles * (num_datasets // len(linestyles) + 1)
            linestyles = linestyles[:num_datasets]
        else:
            linestyles = ['-'] * num_datasets
            
        # Handle line widths - allow both single value and list
        if linewidths is not None:
            linewidths = linewidths * (num_datasets // len(linewidths) + 1)
            linewidths = linewidths[:num_datasets]
        else:
            linewidths = [linewidth] * num_datasets
            
        if legend_labels is not None:
            if len(legend_labels) < num_datasets:
                legend_labels.extend([f'Data {i+1}' for i in range(len(legend_labels), num_datasets)])
        else:
            legend_labels = [f'Data {i+1}' for i in range(num_datasets)]
        
        # Plot each dataset
        for i, dataset in enumerate(data):
            x, y = dataset[:, 0], dataset[:, 1]
            ax.plot(x, y, 
                   color=colors[i], 
                   linestyle=linestyles[i], 
                   linewidth=linewidths[i], 
                   label=legend_labels[i])

    def _set_plot_labels(self, ax, title, xlabel, ylabel, normalized, title_fontsize=24, title_weight='bold', title_font='Arial', xlabel_fontsize=24, xlabel_weight='normal', xlabel_font='Arial', ylabel_fontsize=24, ylabel_weight='normal', ylabel_font='Arial'):
        """Helper method to set plot labels with proper formatting."""
        if title:
            ax.set_title(title, fontsize=title_fontsize, fontweight=title_weight, fontfamily=title_font, pad=20)
            
        # Set x-axis label
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=xlabel_fontsize, fontweight=xlabel_weight, fontfamily=xlabel_font, labelpad=10)
        else:
            ax.set_xlabel('Volume (ml)', fontsize=xlabel_fontsize, fontweight=xlabel_weight, fontfamily=xlabel_font, labelpad=10)
            
        # Set y-axis label
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=ylabel_fontsize, fontweight=ylabel_weight, fontfamily=ylabel_font, labelpad=10)
        else:
            ax.set_ylabel('Normalized Absorbance' if normalized else 'Absorbance (mAU)', 
                        fontsize=ylabel_fontsize, fontweight=ylabel_weight, fontfamily=ylabel_font, labelpad=10)

    def _configure_axes(self, ax, normalized, x_major_ticks, y_major_ticks, x_minor_ticks, 
                       y_minor_ticks, major_tick_length, minor_tick_length, major_tick_width, minor_tick_width, show_box):
        """Configure plot axes settings."""
        # X-axis ticks
        if x_major_ticks:
            ax.xaxis.set_major_locator(MultipleLocator(x_major_ticks))
        if x_minor_ticks:
            ax.xaxis.set_minor_locator(MultipleLocator(x_minor_ticks))
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator())

        # Y-axis ticks
        if normalized:
            if y_major_ticks:
                # Use custom major tick spacing for normalized data
                ax.yaxis.set_major_locator(plt.MultipleLocator(y_major_ticks))
                # Generate tick positions between 0 and 1
                major_ticks = np.arange(0, 1.1, y_major_ticks)
                ax.set_yticks(major_ticks)
                # Format tick labels with consistent decimal places
                ax.set_yticklabels([f'{tick:0.1f}' for tick in major_ticks])
            else:
                # Default to 0.2 intervals if no custom ticks specified
                ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
                major_ticks = np.arange(0, 1.1, 0.2)
                ax.set_yticks(major_ticks)
                ax.set_yticklabels([f'{tick:0.1f}' for tick in major_ticks])
            
            if y_minor_ticks:
                ax.yaxis.set_minor_locator(plt.MultipleLocator(y_minor_ticks))
            else:
                ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
            
            y_padding = 0.05
            ax.set_ylim(-y_padding, 1 + y_padding)
        else:
            if y_major_ticks:
                ax.yaxis.set_major_locator(MultipleLocator(y_major_ticks))
            if y_minor_ticks:
                ax.yaxis.set_minor_locator(MultipleLocator(y_minor_ticks))
            else:
                ax.yaxis.set_minor_locator(AutoMinorLocator())

        # Tick parameters
        ax.tick_params(axis='both', which='major', labelsize=24, length=major_tick_length, 
                      direction='out', top=True, bottom=True, left=True, right=True, width=major_tick_width)
        ax.tick_params(axis='both', which='minor', length=minor_tick_length, 
                      direction='out', top=True, bottom=True, left=True, right=True, width=minor_tick_width)

        # Make tick labels bold
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_weight('bold')

        # Box visibility
        for spine in ax.spines.values():
            spine.set_visible(show_box)
            if show_box:
                spine.set_linewidth(2.0)  # Make the border thicker
        ax.margins(x=0.02, y=0.02)

    def stacked_plot(self, groups, data=None, normalized=False,
                     title='', xlabel='', ylabel='',
                     legend_labels=None, colors=None, linestyles=None,
                     linewidth=2.5, linewidths=None,
                     x_range=None, y_range=None, y_ranges=None,
                     figsize=None, subplot_height=4, dpi=300, output_file=None,
                     hspace=0.0,
                     x_major_ticks=None, y_major_ticks=None,
                     x_minor_ticks=None, y_minor_ticks=None,
                     major_tick_length=10, minor_tick_length=6,
                     major_tick_width=2, minor_tick_width=1,
                     show_box=True, legend_box=False,
                     legend_position=None, legend_fontsize=24,
                     title_fontsize=24, title_weight='bold', title_font='Arial',
                     xlabel_fontsize=24, xlabel_weight='normal', xlabel_font='Arial',
                     ylabel_fontsize=24, ylabel_weight='normal', ylabel_font='Arial',
                     ylabel_per_subplot=None):
        """
        Create a publication-quality stacked Matplotlib plot.

        Each element of `groups` is a list of trace indices to overlay in one
        subplot. Subplots are stacked vertically with a shared x-axis.

        Examples:
            groups=[[0], [1], [2]]   -> 3 subplots, 1 trace each (fully stacked)
            groups=[[0, 1], [2]]     -> 2 subplots: overlay traces 0 & 1, then trace 2 alone
            groups=[[0], [1, 2]]     -> 2 subplots: trace 0 alone, then overlay traces 1 & 2

        Args:
            groups (list of lists of int): Stack structure (see examples above).
            data, normalized, x_range, x_*_ticks, *_tick_*, show_box, etc.:
                Same meaning as in list_plot. Styling args (colors, linestyles,
                linewidths, legend_labels) are indexed by global trace index, so
                trace i always has the same color regardless of which subplot it
                lands in.
            y_ranges (list): Per-subplot y-axis ranges, e.g. [(0, 1), (0, 100)].
                Use None inside the list to fall back to y_range / auto for
                that subplot. Falls back to y_range if not provided.
            subplot_height (float): Height (inches) for each subplot when
                figsize is None.
            hspace (float): Vertical spacing between subplots. 0 = touching
                (shared x-axis tick labels hidden on upper subplots).
            ylabel_per_subplot (bool or None): If True, label every subplot's
                y-axis; if False, only label the middle subplot. None (default)
                auto-picks: False when hspace == 0 (avoids overlapping labels
                on flush-stacked subplots), True otherwise.
        """
        plt.rcParams['font.family'] = 'Arial'

        # Resolve data into a flat list of [N, 2] arrays
        if data is None:
            all_arrays = self.process_data_to_arrays(normalize=normalized)
        elif isinstance(data, SECData):
            all_arrays = data.process_data_to_arrays(normalize=normalized)
        elif isinstance(data, list):
            all_arrays = data
        else:
            all_arrays = [data]

        n_traces = len(all_arrays)

        # Validate groups
        if not groups or not all(isinstance(g, (list, tuple)) and len(g) > 0 for g in groups):
            raise ValueError("groups must be a non-empty list of non-empty lists, e.g. [[0,1],[2]]")
        for grp in groups:
            for idx in grp:
                if not isinstance(idx, (int, np.integer)) or idx < 0 or idx >= n_traces:
                    raise ValueError(f"Group index {idx} out of bounds (have {n_traces} traces)")

        n_subplots = len(groups)

        # Auto-resolve y-label placement: avoid overlapping labels when stacked flush
        if ylabel_per_subplot is None:
            ylabel_per_subplot = hspace > 0 or n_subplots == 1

        # Per-trace style arrays (indexed by global trace index)
        if colors is None:
            colors = list(plt.cm.viridis(np.linspace(0, 1, n_traces)))
        else:
            colors = (list(colors) * (n_traces // len(colors) + 1))[:n_traces]

        if linestyles is None:
            linestyles = ['-'] * n_traces
        else:
            linestyles = (list(linestyles) * (n_traces // len(linestyles) + 1))[:n_traces]

        if linewidths is None:
            linewidths = [linewidth] * n_traces
        else:
            linewidths = (list(linewidths) * (n_traces // len(linewidths) + 1))[:n_traces]

        if legend_labels is None:
            legend_labels = [f'Data {i+1}' for i in range(n_traces)]
        elif len(legend_labels) < n_traces:
            legend_labels = list(legend_labels) + [
                f'Data {i+1}' for i in range(len(legend_labels), n_traces)
            ]

        # Build the figure
        if figsize is None:
            figsize = (12, subplot_height * n_subplots)

        fig, axes = plt.subplots(n_subplots, 1, figsize=figsize, dpi=dpi, sharex=True)
        if n_subplots == 1:
            axes = [axes]

        # Plot each subplot
        for sp_idx, (ax, group) in enumerate(zip(axes, groups)):
            for trace_idx in group:
                arr = all_arrays[trace_idx]
                x, y = arr[:, 0], arr[:, 1]
                ax.plot(x, y,
                        color=colors[trace_idx],
                        linestyle=linestyles[trace_idx],
                        linewidth=linewidths[trace_idx],
                        label=legend_labels[trace_idx])

            # Axis configuration (ticks, box, normalized y-range default)
            self._configure_axes(ax, normalized, x_major_ticks, y_major_ticks,
                                 x_minor_ticks, y_minor_ticks,
                                 major_tick_length, minor_tick_length,
                                 major_tick_width, minor_tick_width, show_box)

            # Per-subplot y-range overrides
            sp_y_range = None
            if y_ranges is not None and sp_idx < len(y_ranges):
                sp_y_range = y_ranges[sp_idx]
            if sp_y_range is None and y_range is not None:
                sp_y_range = y_range
            if sp_y_range is not None:
                ax.set_ylim(sp_y_range)

            # Y-axis label
            if ylabel_per_subplot or sp_idx == n_subplots // 2:
                ylabel_text = ylabel or ('Normalized Absorbance' if normalized else 'Absorbance (mAU)')
                ax.set_ylabel(ylabel_text, fontsize=ylabel_fontsize,
                              fontweight=ylabel_weight, fontfamily=ylabel_font, labelpad=10)

            # Per-subplot legend
            if legend_position:
                leg = ax.legend(fontsize=legend_fontsize, frameon=legend_box,
                                bbox_to_anchor=legend_position)
            else:
                leg = ax.legend(fontsize=legend_fontsize, frameon=legend_box)
            if leg:
                for text in leg.get_texts():
                    text.set_weight('bold')

        # Shared x-axis range and bottom label
        if x_range:
            for ax in axes:
                ax.set_xlim(x_range)

        xlabel_text = xlabel or 'Volume (ml)'
        axes[-1].set_xlabel(xlabel_text, fontsize=xlabel_fontsize,
                            fontweight=xlabel_weight, fontfamily=xlabel_font, labelpad=10)

        # Title on top subplot
        if title:
            axes[0].set_title(title, fontsize=title_fontsize, fontweight=title_weight,
                              fontfamily=title_font, pad=20)

        plt.subplots_adjust(hspace=hspace)

        # When stacked flush, hide x tick labels on every subplot but the last
        if hspace == 0:
            for ax in axes[:-1]:
                plt.setp(ax.get_xticklabels(), visible=False)

        if hspace > 0:
            plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
            print(f"Figure saved as {output_file}")
        else:
            plt.show()

        return fig, axes

    def integrate_peaks(self, regions, baseline_corrected=True):
        """
        Calculate area under the curve for specified regions (e.g., monomer, oligomer peaks).
        
        Args:
            regions (dict): Dictionary of peak regions with format:
                          {'peak_name': (start_vol, end_vol)}
                          e.g., {'monomer': (12, 15), 'oligomer': (8, 11)}
            baseline_corrected (bool): Whether to use baseline-corrected data
            
        Returns:
            list of dicts: List containing AUC values for each dataset and region
        """
        results = []
        
        for i, df in enumerate(self.data):
            peak_areas = {'dataset': i}
            
            # Choose which absorbance column to use
            y_col = 'Adjusted_Absorbance_mAU' if baseline_corrected and 'Adjusted_Absorbance_mAU' in df.columns else 'Absorbance_mAU'
            
            for peak_name, (start_vol, end_vol) in regions.items():
                # Get data within peak region
                mask = (df['Volume_ml'] >= start_vol) & (df['Volume_ml'] <= end_vol)
                peak_data = df[mask]
                
                if not peak_data.empty:
                    # Calculate area using trapezoidal integration
                    area = simpson(peak_data[y_col], peak_data['Volume_ml'])
                    peak_areas[f'{peak_name}_area'] = area
                    
                    # Calculate peak height and position
                    max_idx = peak_data[y_col].idxmax()
                    peak_areas[f'{peak_name}_height'] = peak_data.loc[max_idx, y_col]
                    peak_areas[f'{peak_name}_position'] = peak_data.loc[max_idx, 'Volume_ml']
                else:
                    peak_areas[f'{peak_name}_area'] = 0
                    peak_areas[f'{peak_name}_height'] = 0
                    peak_areas[f'{peak_name}_position'] = None
            
            results.append(peak_areas)
        
        return results

    def calculate_monomer_ratio(self, monomer_region, total_protein_region=None, baseline_corrected=True):
        """
        Calculate the ratio of monomer to total protein.
        
        Args:
            monomer_region (tuple): (start_vol, end_vol) for monomer peak
            total_protein_region (tuple): Optional (start_vol, end_vol) for total protein.
                                        If None, uses full range
            baseline_corrected (bool): Whether to use baseline-corrected data
            
        Returns:
            list: Monomer ratios for each dataset
        """
        regions = {'monomer': monomer_region}
        if total_protein_region:
            regions['total'] = total_protein_region
            
        results = []
        areas = self.integrate_peaks(regions, baseline_corrected)
        
        for i, area_data in enumerate(areas):
            if total_protein_region:
                total_area = area_data['total_area']
            else:
                # Use full range if total_protein_region not specified
                df = self.data[i]
                y_col = 'Adjusted_Absorbance_mAU' if baseline_corrected and 'Adjusted_Absorbance_mAU' in df.columns else 'Absorbance_mAU'
                total_area = simpson(df[y_col], df['Volume_ml'])
            
            monomer_ratio = area_data['monomer_area'] / total_area if total_area > 0 else 0
            results.append({
                'dataset': i,
                'monomer_area': area_data['monomer_area'],
                'total_area': total_area,
                'monomer_ratio': monomer_ratio
            })
            
        return results
