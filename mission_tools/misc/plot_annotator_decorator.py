from datetime import datetime
import matplotlib.pyplot as plt
import pathlib
import inspect

def annotate_plot(func):
    def wrapper_annotate_plot(*args, **kwargs):
        func(*args, **kwargs)
        plot_date_time = datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S")

        # Get the filename assisiated with func.
        frame_info = inspect.stack()
        filepath = frame_info[1][1]
        file_name = pathlib.Path(filepath).name

        plt.subplots_adjust(bottom=0.2)
        plt.text(10, 10, 
                f'Generated at {plot_date_time} by {func.__name__}() in {file_name}', 
                transform=None)
    return wrapper_annotate_plot