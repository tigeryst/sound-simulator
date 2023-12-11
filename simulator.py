# %% [markdown]
# # Simulator

# ## Background
# This computer simulation aims to numerically simulate the behaviour of acoustic
# systems involving arrays of transducers (emitters and receivers) positioned on a plane.
#
# Each emitter acts as a spherical point source emitting sinusoidal waves of defined
# frequency and amplitude for a number of cycles. Sound propagates at the defined speed
# of sound and frequency. Its amplitude is attenuated by the inverse of the distance
# between the emitter and the receiver.
#
# Given the readings from an array of receivers, the simulator can be used to reconstruct
# the acoustic sources within a defined region of interest using a beamforming algorithm.


# ## Contents
# This file implements the following classes for the acoustic simulator:
# 1. Transducer
# 1. Receiver
# 1. Emitter
# 1. SoundSimulator
# 1. BeamFormer

# %% [markdown]
# ## Implementation

# ### Imports
# %%
import math


# %% [markdown]
# ### Emitter and Receiver
# %%
class Transducer:
    """Device to convert energy from one form to the other.

    Attributes:
        x (float): x-coordinate in metres of the transducer.
        y (float): y-coordinate in metres of the transducer.
        t_array (List[float]): Array of time points in seconds where the signal is sampled.
        signal (List[float]): Array of displacement at each time point.
    """

    def __init__(self, x: float, y: float, t_array: list):
        """Initialize the transducer with its position and time array.

        Args:
            x (float): x-coordinate in metres of the transducer.
            y (float): y-coordinate in metres of the transducer.
            t_array (list): Array of time points in seconds where the signal is sampled.
        """
        self.x = x
        self.y = y
        self.t_array = t_array
        self.signal = len(self.t_array) * [0.0]


# %%
class Receiver(Transducer):
    """Transducer that converts acoustic energy into electrical energy.

    Attributes:
        x (float): x-coordinate in metres of the receiver.
        y (float): y-coordinate in metres of the receiver.
        t_array (List[float]): Array of time points in seconds where the signal is sampled.
        signal (List[float]): Array of displacement at each time point.
    """

    def __init__(self, x: float, y: float, t_array: list):
        """Initialise the receiver with its position and signal array.

        Args:
            x (float): x-coordinate in metres of the receiver.
            y (float): y-coordinate in metres of the receiver.
            t_array (list): Array of time points in seconds where the signal is sampled.
        """
        super().__init__(x, y, t_array)


class Emitter(Transducer):
    """Transducer that converts electrical energy into acoustic energy.

    Attributes:
        x (float): x-coordinate in metres of the emitter.
        y (float): y-coordinate in metres of the emitter.
        t_array (List[float]): Array of time points in seconds where the signal is sampled.
        signal (List[float]): Array of displacement at each time point.

    Methods:
        generate_signal: Generate a sinusoidal signal and store it in the signal array.
    """

    def __init__(self, x: float, y: float, t_array: list):
        """Initialise the emitter with its position and signal array.

        Args:
            x (float): x-coordinate in metres of the emitter.
            y (float): y-coordinate in metres of the emitter.
            t_array (list): Array of time points in seconds where the signal is sampled.
        """
        super().__init__(x, y, t_array)

    def generate_signal(self, f_c: float, n_cycles: int, amplitude: float):
        """Generate a sinusoidal signal and store it in the signal array.

        Args:
            f_c (float): Frequency of the sinusoidal signal in Hz.
            n_cycles (int): Number of cycles of the sinusoidal signal.
            amplitude (float): Amplitude of the sinusoidal signal.

        Returns:
            signal (List[float]): Array of displacement at each time point.
        """
        # Theory:
        # Define a sinusoidal function within the given domain assuming zero displacement at t = t_min
        # s = A * sin[2π * f_c * (t - t_min)] for t_min <= t <= t_max
        # s = 0 otherwise
        # t_min = t_array[0]
        # t_max = n_cycles*T where the time period T = 1/f_c

        time_period = 1 / f_c
        t_min = self.t_array[0]
        t_max = n_cycles * time_period
        self.signal = [
            amplitude * math.sin(2 * math.pi * f_c * (t - t_min))
            if t_min <= t <= t_max
            else 0.0
            for t in self.t_array
        ]
        return self.signal


class SoundSimulator:
    """Simulate the propagation of sound in a medium.

    Attributes:
        emitters (List[Emitter]): Array of emitters.
        receivers (List[Receiver]): Array of receivers.
        t_array (List[float]): Array of time points in seconds where the signal is sampled.
        sos (float): Speed of sound in the medium in metres per second.

    Methods:
        run: Simulate the propagation of sound from emitters to receivers.
    """

    def __init__(
        self,
        emitters: list = None,
        receivers: list = None,
        t_array: list = None,
        sos: float = 1500.0,
    ):
        """Initialise the simulator with the given parameters.

        Args:
            emitters (List[Emitter], optional): Array of emitters. Defaults to None.
            receivers (List[Receiver], optional): Array of receivers. Defaults to None.
            t_array (List[float], optional): Array of time points in seconds where the signal is sampled. Defaults to None.
            sos (float, optional): Speed of sound in the medium in metres per second. Defaults to 1500.0.
        """
        self.emitters = emitters if emitters else []
        self.receivers = receivers if receivers else []
        self.t_array = t_array if t_array else []
        self.sos = sos

    def run(self):
        """Simulate the propagation of sound from emitters to receivers.

        Returns:
            simulated_receivers (List[Receiver]): Array of receivers with their respective signal arrays updated to reflect the result of the superposition of sound waves from the array of emitters.
        """
        simulated_receivers = []
        for receiver in self.receivers:
            receiver_signal = [0.0] * len(self.t_array)

            for emitter in self.emitters:
                distance = math.sqrt(
                    (receiver.x - emitter.x) ** 2 + (receiver.y - emitter.y) ** 2
                )
                time_delay = distance / self.sos

                for i, t in enumerate(self.t_array):
                    if t >= time_delay:
                        attenuation = 1 / distance
                        receiver_signal[i] += (
                            emitter.signal[
                                i
                                - int(time_delay / (self.t_array[1] - self.t_array[0]))
                            ]
                            * attenuation
                        )

            simulated_receiver = Receiver(receiver.x, receiver.y, self.t_array)
            simulated_receiver.signal = receiver_signal
            simulated_receivers.append(simulated_receiver)

        self.receivers = simulated_receivers
        return simulated_receivers
