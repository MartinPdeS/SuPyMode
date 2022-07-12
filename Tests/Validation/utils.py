from mayavi              import mlab
import matplotlib._pylab_helpers
import matplotlib.pyplot as plt
from pyface.api           import GUI

def CloseMatplotlib():
    """Close matplotlib scene."""
    figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    if len(figures) >= 1:
        plt.close('all')


def CloseMlab():
    """Close mayavi scene."""
    engine = mlab.get_engine()
    scene = engine.current_scene
    if scene is not None:
        mlab.close()

def Close():
    """Close all scene."""
    CloseMatplotlib()

    CloseMlab()


def TestPlot(func):
    def Test(func):
        def wrapper(*args):
            GUI.invoke_after(5000, Close)
            func(*args)

        return wrapper
    return Test



def TestFactory(Name):
    def Test(func):
        def wrapper(*args):
            try:
                func(*args)
                print(f'       {func.__name__}: {Name}.'.ljust(100) + '-'*15 + '->  PASSED')
            except Exception as error:
                print(f'       {func.__name__}: {Name}.'.ljust(100) + '-'*15 + f'->  FAILED\n{error}')

        return wrapper
    return Test


class BaseStepTest():
    def _steps(self):
        for name in dir(self):
            if name.startswith("step"):
                yield name, getattr(self, name)

    def test_steps(self):
        for name, step in self._steps():
            try:
                step()
            except Exception as e:
                self.fail("{} failed ({}: {})".format(step, type(e), e))
