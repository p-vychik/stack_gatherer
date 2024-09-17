from multiprocessing import Pipe
from juliacall import Main as jl


def compute_avg_speed(child_conn):
    while True:
        if child_conn.poll():
            message = child_conn.recv()
            if isinstance(message, str):
                jl.include(message)
                continue
            m_1, m_2 = message
            avg_speed = jl.fn(m_1, m_2)
            if isinstance(avg_speed[-1], float) or isinstance(avg_speed[-1], int):
                avg_speed = round(avg_speed[-1], 3)
                child_conn.send(avg_speed)
            else:
                child_conn.send(None)
