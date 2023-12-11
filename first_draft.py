import math


class Transducer:
    def __init__(self, x, y, t_array):
        self.x = x
        self.y = y
        self.t_array = t_array
        self.signal = len(self.t_array) * [0]


class Receiver(Transducer):
    def __init__(self, x, y, t_array):
        super().__init__(x, y, t_array)


class Emitter(Transducer):
    def __init__(self, x, y, t_array):
        super().__init__(x, y, t_array)

    def generate_signal(self, f_c, n_cycles, amplitude):
        # Calculate the time step from the first two elements in t_array
        dt = self.t_array[1] - self.t_array[0]

        # Generate a sinusoidal signal
        signal_time = [i * dt for i in range(int(n_cycles / f_c / dt))]
        sinusoidal_signal = [
            amplitude * math.sin(2 * math.pi * f_c * t) for t in signal_time
        ]

        # Find the index where the sinusoidal signal should start in self.t_array
        start_index = int(self.t_array[0] / dt)

        # Initialize self.signal with zeros
        self.signal = len(self.t_array) * [0]

        # Copy the generated sinusoidal signal into the appropriate part of self.signal
        for i in range(len(sinusoidal_signal)):
            if start_index + i < len(self.signal):
                self.signal[start_index + i] = sinusoidal_signal[i]

        # Return the generated signal
        return self.signal


class SoundSimulator:
    def __init__(self, emitters=None, receivers=None, t_array=None, sos=1500.0):
        self.emitters = emitters if emitters is not None else []
        self.receivers = receivers if receivers is not None else []
        self.t_array = t_array if t_array is not None else []
        self.sos = sos

    def run(self):
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


import matplotlib.pyplot as pyplot


# from simulator import *
def main():
    # ---------------------------------------------------------------
    # setup the 'clock' with which the simulation will be run
    t_delta = 0.1e-6  # time step size [seconds]
    t_N = 800  # number of time points
    t_array = [
        t_i * t_delta for t_i in range(t_N)
    ]  # the master clock, array of time [seconds]
    sos = 1500.0  # speed of sound [metres/second]
    # ---------------------------------------------------------------
    # setup receivers and emitters
    # arrange a linear order of receivers from -0.04 m to 0.04 m
    receiver_positions = [
        [receiver_x / 1000.0, 0.04] for receiver_x in range(-40, 41, 2)
    ]
    emitters = []
    receivers = []
    emitters.append(Emitter(0, 0, t_array))
    signal = emitters[0].generate_signal(1e6, 5, 1)
    fig, axs = pyplot.subplots(1)
    axs.plot(t_array, signal)
    axs.set_xlabel("time (seconds)")
    axs.set_ylabel("amplitude")
    pyplot.show()

    for pos in receiver_positions:
        receivers.append(Receiver(pos[0], pos[1], t_array))
    # ---------------------------------------------------------------
    # run the simulation
    simulation = SoundSimulator(emitters, receivers, t_array, sos)
    receivers = simulation.run()
    img = []
    for i in range(len(receivers)):
        img.append(receivers[i].signal)
    fig, axs = pyplot.subplots(1)
    imgplot = axs.imshow(img)
    axs.set_aspect(20)
    axs.set_xlabel("time point")
    axs.set_ylabel("receiver number")
    fig.colorbar(imgplot)
    pyplot.show()


if __name__ == "__main__":
    main()
