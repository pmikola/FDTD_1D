import math as m


def source(field, current_time, delay, source_width, denom):
    if field == "Hx":
        return m.exp(-(((current_time - delay) * (current_time - delay)) / source_width / denom))
    if field == "Ey":
        return m.exp(-(((current_time - delay) * (current_time - delay))/source_width / denom))
        # return m.exp(-0.5 * (m.pow((current_time - delay) / source_width, 2)))


def source_sin(frequency, source_width, dt):
    return m.sin(2 * m.pi * frequency * dt * source_width)
