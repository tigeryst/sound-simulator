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
        """Initialise the transducer with its position and time array.

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
        run: Simulate the propagation of sound from emitters to receivers and update the signal array of the receivers.
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
        """Simulate the propagation of sound from emitters to receivers and update the signal array of the receivers.

        Receivers receive the superposition of sound waves from all emitters.
        The amplitude of the sound wave is attenuated by the inverse of the distance between the emitter and the receiver.
        The signal is delayed by the time it takes for the sound wave to travel from the emitter to the receiver.
        The signal at each receiver is stored in the signal array of the receiver.

        Returns:
            simulated_receivers (List[Receiver]): Array of receivers with their respective signal arrays updated to
            reflect the result of the superposition of sound waves from the array of emitters.
        """
        # Theory for attenuation:
        # Amplitude is proportional to the inverse of the distance between the emitter and the receiver.
        # Assume that the amplitude of the signal at the emitter, A_0, is measured at a reference distance, d_0.
        # Assume the far-field approximation where the distance between the emitter and the receiver is much
        # greater than the size of the emitter and greater than d_0.
        # A = (d_0 / d) * A_0 where d is the distance between the emitter and the receiver.
        # Assume d_0 = 1m for this exercise.

        # Theory for time delay:
        # The time delay is given by the distance between the emitter and the receiver divided by the speed of sound.

        # Assume that t_array for all transducers are the same as that defined in the SoundSimulator object.
        n = len(self.t_array)  # number of time points/signal samples
        dt = self.t_array[1] - self.t_array[0]  # time step size
        for receiver in self.receivers:
            receiver.signal = n * [0.0]
            # Shift the signal from each emitter by the time it takes to travel from the emitter to the receiver.
            for emitter in self.emitters:
                # Calculate the time it takes for the sound wave to travel from the emitter to the receiver.
                distance = math.dist([receiver.x, receiver.y], [emitter.x, emitter.y])
                delay_time = distance / self.sos
                # ? Confirm with supervisor whether the time step can be assumed to be constant.
                # * Assume that the time step size is constant as would normally be when sampling signals.
                # * Assume that the signals are properly sampled such that the time step size is much
                # * smaller than the time period of the signal.
                # Find the number of time steps to shift the signal by.
                delay_index = int(delay_time / dt)
                for i in range(n):
                    if i >= delay_index:
                        # Add the time-delayed, attenuated signal from the emitter to the receiver.
                        receiver.signal[i] += emitter.signal[i - delay_index] / distance
        # Return the receivers with their signal arrays updated.
        return self.receivers


class BeamFormer:
    """Reconstruct the acoustic sources within a region of interest.

    Attributes:
        receivers (List[Receiver]): Array of receivers.
        x_array (List[float]): Array of x-coordinates of the region of interest.
        y_array (List[float]): Array of y-coordinates of the region of interest.
        t_array (List[float]): Array of time points in seconds where the signal is sampled.
        sos (float): Speed of sound in the medium in metres per second.
        field (List[List[List[float]]]): Beamformed signal at each point in the region of interest and each point in time.

    Methods:
        generate_field: Reconstruct the acoustic sources within the region of interest.
    """

    def __init__(
        self, receivers: list, x_array: list, y_array: list, t_array: list, sos: float
    ):
        """Initialise the beamformer with the given parameters.

        Args:
            receivers (List[Receiver]): Array of receivers.
            x_array (List[float]): Array of x-coordinates of the region of interest.
            y_array (List[float]): Array of y-coordinates of the region of interest.
            t_array (List[float]): Array of time points in seconds where the signal is sampled.
            sos (float): Speed of sound in the medium in metres per second.
        """
        self.receivers = receivers
        self.x_array = x_array
        self.y_array = y_array
        self.t_array = t_array
        self.sos = sos
        self.field = [
            [[0.0 for _ in range(len(self.t_array))] for _ in range(len(self.x_array))]
            for _ in range(len(self.y_array))
        ]

    def generate_field(self):
        """Reconstruct the acoustic sources within the region of interest using the delay-and-sum beamforming method.

        The reconstructed field is stored in the field array.
        """
        # Theory:
        # The signal at each point in the region of interest is the weighted superposition of the signals detected at each receiver.
        # Signal from the point arrrives at each receiver at different times. The delay-and-sum method involves delaying the signal
        # from each receiver to account for this time difference and superposing the signals such that they are in phase.
        # Components of the detected signal originating from sources other than the point of interest will be out of phase and suppressed.

        # We also reverse the effects of attenuation by multiplying the signal from each receiver by the inverse of the attenuation factor

        # q(r, t) = 1/N * sum_{i=1}^{N} [ (d_i/d_0) * p(r_i, t - tau_i) ] where N is the number of receivers, r is the point of interest,
        # r_i is the position of the i-th receiver, d_i is the distance between the i-th receiver and the point of interest, p is the
        # signal detected at the i-th receiver, and tau_i is the time delay between the signal from the i-th receiver and the point of
        # interest. The time delay is given by tau_i = d_i / c where c is the speed of sound. Assume d_0 = 1m.

        n = len(self.receivers)  # number of receivers

        # Iterate through each point of interest
        for i, x in enumerate(self.x_array):
            for j, y in enumerate(self.y_array):
                # Iterate through each receiver:
                for receiver in self.receivers:
                    # Calculate the distance between the receiver and the point of interest.
                    distance = math.dist([x, y], [receiver.x, receiver.y])
                    # Calculate the time delay between the signal from the receiver and the point of interest.
                    delay_time = distance / self.sos

                    # Find the number of time steps to shift the signal by.
                    delay_index = int(delay_time / (self.t_array[1] - self.t_array[0]))
                    # Add the time-delayed, attenuated signal from the receiver to the point of interest.
                    for k in range(len(self.t_array)):
                        if k + delay_index < len(self.t_array):
                            self.field[j][i][k] += (
                                distance * receiver.signal[k + delay_index] / n
                            )
                    # if i == 0 and j == 0:
                    #     print(delay_index)
                    # for k in range(len(self.t_array)):
                    #     if k >= delay_index:
                    #         self.field[j][i][k] += (
                    #             distance * receiver.signal[k - delay_index] / n
                    #         )
