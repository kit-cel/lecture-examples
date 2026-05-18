import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import sounddevice as sd


# ----- Quantisierer-Definition ------------------------------------------
def midtread_quantizer(x,w,x_max):
    Delta_x = x_max / (2**(w-1));
    xh = np.sign(x)*Delta_x*np.floor(np.abs(x)/Delta_x+0.5);   #Select class from +/- 2**w-1/2
    xh = np.clip(xh, -1*(2**(w-1)-1)*Delta_x, (2**(w-1)-1)*Delta_x)
    return xh


FS = 44100                                 # Abtastrate [Hz]
DURATION = 1.0                             # Länge des Test‑Signals [s]

F_MIN = 50
F_MAX = 4000

# --------------------------------------------------------------
# Hilfsfunktion: erzeugt Signal für gegebene Frequenz
# --------------------------------------------------------------
#def gen_signal(freq):
#    """
#    freq      – Frequenz in Hz
#    """
#    t = np.linspace(0, DURATION, int(FS * DURATION), endpoint=False)
#    return 0.4 * np.sin(2 * np.pi * freq * t)
#   

# --------------------------------------------------------------
# Klasse, die Figure, Plot, Slider und Play‑Button kapselt
# --------------------------------------------------------------
class AudioFigure:
    def __init__(self, init_freq, init_resolution, title):
        self.freq = init_freq
        self.w = init_resolution
        self.sr = FS
        self.fs = 2*init_freq

        # ---- Plot und Daten ------------------------------------
        self.fig, self.ax = plt.subplots(2,2)
        self.fig.subplots_adjust(bottom=0.30)   # Platz für Slider + Button
        self.gen_signal()
        self.line1, = self.ax[0,0].plot(
            self.t[0:1000],
            self.samples[0:1000],
            color='tab:blue'
        )
        self.ax[0,0].set_xlabel('Zeit [s]')
        self.ax[0,0].set_ylabel('Amplitude')
        self.ax[0,0].set_title(title)
        self.ax[0,0].set_ylim(1.1*np.min(self.samples),1.1*np.max(self.samples))

        self.line2, = self.ax[0,1].plot(
            self.t[0:1000],
            self.samples_q[0:1000],
            color='tab:blue'
        )
        self.ax[0,1].set_xlabel('Zeit [s]')
        self.ax[0,1].set_ylabel('Amplitude')
        self.ax[0,1].set_title(title+" quantisiert")
        self.ax[0,1].set_ylim(1.1*np.min(self.samples),1.1*np.max(self.samples))

        #Spektrum
        self.line3, = self.ax[1,0].plot(
            np.linspace(-self.sr/2,self.sr/2,self.samples.size),
            10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)),
            color='tab:blue'
        )
        self.ax[1,0].set_xlabel('Freq [Hz]')
        self.ax[1,0].set_ylabel('Betragsquadrat')
        self.ax[1,0].set_title(title+" Spektrum")
        self.ax[1,0].set_xlim(-1.1*F_MAX,1.1*F_MAX)
        self.ax[1,0].set_ylim(-100,0)

        #Spektrum
        self.line4, = self.ax[1,1].plot(
            np.linspace(-self.sr/2,self.sr/2,self.samples.size),
            10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)),
            color='tab:blue'
        )
        self.ax[1,1].set_xlabel('Freq [Hz]')
        self.ax[1,1].set_ylabel('Betragsquadrat')
        self.ax[1,1].set_title(title+" Spektrum")
        self.ax[1,1].set_ylim(-100,0)
        

        # ---- Play‑Buttons ---------------------------------------
        btn_ax1 = self.fig.add_axes([0.2, 0.10, 0.2, 0.1])     # [l, b, w, h] in Fig‑Koord.
        self.play_btn1 = Button(
            btn_ax1, '▶ Play',
            color='lightgoldenrodyellow',
            hovercolor='0.975'
        )
        self.play_btn1.on_clicked(self._play1)

        btn_ax2 = self.fig.add_axes([0.6, 0.10, 0.2, 0.1])     # [l, b, w, h] in Fig‑Koord.
        self.play_btn2 = Button(
            btn_ax2, '▶ Play',
            color='lightgoldenrodyellow',
            hovercolor='0.975'
        )
        self.play_btn2.on_clicked(self._play2)

        # ---- Frequenz‑Slider ------------------------------------
        slider_f_ax = self.fig.add_axes([0.15, 0.04, 0.70, 0.05]) # [l, b, w, h] in Fig‑Koord.
        self.slider_f = Slider(
            slider_f_ax,
            label='Frequenz [Hz]',
            valmin=F_MIN,
            valmax=F_MAX,
            valinit=self.freq,
            valstep=1,
            color='tab:orange'
        )
        self.slider_f.on_changed(self._update_signal_f)

        # ---- Resolution‑Slider ------------------------------------
        slider_w_ax = self.fig.add_axes([0.15, 0.00, 0.70, 0.05]) # [l, b, w, h] in Fig‑Koord.
        self.slider_w = Slider(
            slider_w_ax,
            label='Resolution [bit]',
            valmin=2,
            valmax=16,
            valinit=self.w,
            valstep=1,
            color='tab:orange'
        )
        self.slider_w.on_changed(self._update_signal_w)

    # ------------------------------------------------------------------
    # Callback: Play‑Button
    # ------------------------------------------------------------------
    def _play1(self, event):
        sd.play(self.samples.astype(np.float32), self.sr, blocking=False)
    def _play2(self, event):
        sd.play(self.samples_q.astype(np.float32), self.sr, blocking=False)

    # ------------------------------------------------------------------
    # Callback: Slider – Signal neu berechnen und Plot aktualisieren
    # ------------------------------------------------------------------
    def _update_signal_f(self, val):
        self.freq = float(val)
        self.gen_signal()
        # Plot‑Daten ersetzen, ohne neue Figure zu erzeugen
        self.line1.set_ydata(self.samples[0:1000])
        self.line2.set_ydata(self.samples_q[0:1000])
        self.line3.set_ydata(10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)))
        self.line4.set_ydata(10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)))
        self.fig.canvas.draw_idle()

    def _update_signal_w(self, val):
        self.w = float(val)
        self.gen_signal()
        # Plot‑Daten ersetzen, ohne neue Figure zu erzeugen
        self.line1.set_ydata(self.samples[0:1000])
        self.line2.set_ydata(self.samples_q[0:1000])
        self.line3.set_ydata(10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples))+1e-30)**2)))
        self.line4.set_ydata(10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)-10*np.log10(np.max(np.abs(np.fft.fftshift(np.fft.fft(self.samples_q))+1e-30)**2)))
        self.fig.canvas.draw_idle()
    # --------------------------------------------------------------
    # Hilfsfunktion: erzeugt Signal für gegebene Frequenz
    # --------------------------------------------------------------
    def gen_signal(self):
        self.t = np.linspace(0, DURATION, int(FS * DURATION), endpoint=False)    #"Time continious"
        
        self.samples   = 0.4 * np.cos(2 * np.pi * self.freq * self.t)
        DS = np.floor(self.sr/self.fs).astype(int)
        self.samples_q = midtread_quantizer(self.samples,self.w,np.max(self.samples))

    # ------------------------------------------------------------------
    # Öffnen des Fensters (nicht blockierend)
    # ------------------------------------------------------------------
    def show(self):
        self.fig.show()



########################
# main function
########################
def main():

    # --------------------------------------------------------------
    # 2 Fenster erzeugen – jeweils ein eigener Slider
    # --------------------------------------------------------------
    fig_sine = AudioFigure(init_freq=440, init_resolution=4,  title='Sinus')
    #fig_square = AudioFigure(init_freq=660, wave_type='square',
                            #title='Rechteck‑Tone (verstellbare Frequenz)')

    fig_sine.show()
    #fig_square.show()

    # --------------------------------------------------------------
    # Haupt‑Event‑Loop starten (öffnet beide Fenster)
    # --------------------------------------------------------------
    plt.show()


########################
# make it executable
########################
if __name__ == "__main__":
    main()
