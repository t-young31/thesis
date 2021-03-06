

data = [
    [r'\bfseries{R1}',  369.15, 93.3, -44.0, 108.8],
    [r'\bfseries{R2}',  369.15, 77.8, -69.0, 102.2],
    [r'\bfseries{R3}',  322.35, 89.9, -24.5, 98.1 ],
    [r'\bfseries{R4}',  303.15, 66.6, -51.0, 82.0 ],
    [r'\bfseries{R5}',  303.15, 69.1, -28.4, 77.7 ],
    [r'\bfseries{R6}',  303.15, 68.7, -67.3, 89.1 ],
    [r'\bfseries{R7}',  303.15, 77.4, -27.2, 85.6 ],
    [r'\bfseries{R8}',  252.15, 49.8, -51.0, 62.7 ],
    [r'\bfseries{R9}',  303.15, 77.8, -32.6, 87.6 ],
    [r'\bfseries{R10}', 303.15, 84.4, -18.4, 90.0 ]
    ]


if __name__ == '__main__':

    for line in data:
        print(line[0],
              f'{round(line[1], 1):.1f}',
              f'{round(line[2] / 4.184, 1):.1f}',
              f'{round(line[1] * line[3] / (1000 * 4.184), 1):.1f}',  # T∆S
              f'{round(line[4] / 4.184, 1):.1f}',
              sep=' & ')
