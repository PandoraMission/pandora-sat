"""Utilities for plotting, animating, and saving gifs."""
# Standard library
from typing import Optional

# Third-party
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML
from matplotlib import animation as mplanimation
from matplotlib.animation import FFMpegWriter
from PIL import Image

__all__ = ["animate", "save_gif", "save_mp4"]


def _to_matplotlib_animation(
    data,
    instance_name: str = "Frame",
    step: int = None,
    interval: int = 200,
    position: Optional = None,
    figsize: Optional = None,
    **plot_kwargs,
):
    """Lifted from Lightkurve"""

    cmap = plot_kwargs.pop("cmap", "Greys")
    if step is None:
        step = len(data) // 50
        if step < 1:
            step = 1

    if figsize is None:
        x, y = 1.0, 1 * data.shape[1] / data.shape[2]
        if np.max([x, y]) < 5.0:
            c = 5.0 / np.max([x, y])
            x *= c
            y *= c
        x, y = int(x), int(y)
        figsize = (x, y)

    _, ax = plt.subplots(figsize=figsize)
    ax.imshow(data[0], cmap=cmap, **plot_kwargs)
    ax.set(xticks=[], yticks=[])
    if position is not None:
        plt.gca().set_position(position)

    def init():
        return ax.images

    def animate(i):
        frame = i * step
        ax.images[0].set_data(data[frame])
        ax.set_title(f"{instance_name} {frame}")
        return ax.images

    plt.close(ax.figure)
    frames = len(data) // step
    anim = mplanimation.FuncAnimation(
        ax.figure,
        animate,
        init_func=init,
        frames=frames,
        interval=interval,
        blit=True,
    )
    return anim


def animate(data, step: int = None, interval: int = 200, **plot_kwargs):
    """Animate a 3D data set.

    Parameters:
    -----------
    data: np.ndarray
        Data to animate. Has shape (ntime, nrow, ncol), will animate through each frame in ntime.
    step: int
        Which frames should be output. If step is 10, every 10th frame will be show. If step is None,
        step will be calculated such that 50 frames are shown.
    interval:
        Interval between frames in ms.
    """
    return HTML(
        _to_matplotlib_animation(
            data,
            step=step,
            interval=interval,
            position=[0, 0, 1, 1],
            **plot_kwargs,
        ).to_jshtml()
    )


def save_mp4(
    data,
    outfile="out.mp4",
    step: int = None,
    interval: int = 200,
    axis_title="Frame",
    dpi=100,
    **plot_kwargs,
):
    """Create an mp4 from a 3D dataset.

    Parameters:
    -----------
    data: np.ndarray
        Data to animate. Has shape (ntime, nrow, ncol), will animate through each frame in ntime.
    outfile: str
        Path to output file
    step: int
        Which frames should be output. If step is 10, every 10th frame will be show. If step is None,
        step will be calculated such that 50 frames are shown.
    interval:
        Interval between frames in ms.
    axis_title: str,
        Label applied to axis. Default is "Frame".
    dpi: int
        Dots per inch, sets output resolution.
    """
    anim = _to_matplotlib_animation(
        data,
        step=step,
        interval=interval,
        position=[0, 0, 1, 1],
        instance_name=axis_title,
        **plot_kwargs,
    )
    anim.save(
        outfile,
        writer=FFMpegWriter(fps=1000 / interval, bitrate=5000),
        dpi=dpi,
    )


def save_gif(
    data,
    outfile="out.gif",
    step: int = 1,
    interval: int = 50,
    scale: int = 8,
    vmin=-50,
    vmax=50,
):
    """Create a gif from a 3D dataset.

    Parameters:
    -----------
    data: np.ndarray
        Data to animate. Has shape (ntime, nrow, ncol), will animate through each frame in ntime.
    outfile: str
        Path to output file
    step: int
        Which frames should be output. If step is 10, every 10th frame will be show. If step is None,
        step will be calculated such that 50 frames are shown.
    interval:
        Interval between frames in ms.
    scale: int
        gifs are shown at pixel resolution. Use scale to increase the physical size of the gif.
        Setting scale to 8 means that each pixel in `data` is represented by 8x8 pixels on your screen.
    vmin: float
        Minimum color scale value
    vmax: float
        Maximum color scale value
    """
    imgs = (
        (
            np.min(
                [np.max([data, data**0 + vmin], axis=0), data**0 + vmax],
                axis=0,
            )
            - vmin
        )
        * 255
        / (vmax - vmin)
    )
    imgs = imgs.astype(np.uint16)[::step]
    imgs = np.repeat(np.repeat(imgs, scale, axis=1), scale, axis=2)
    imgs = [Image.fromarray(img) for img in imgs]
    imgs[0].save(
        outfile,
        save_all=True,
        append_images=imgs[1:],
        duration=interval,
        loop=0,
    )
