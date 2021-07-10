import gaptrain as gt

if __name__ == '__main__':

    traj = gt.Trajectory(filename='')
    traj.load('GAP_rPBE0-D3_nvt_300K_raw.xyz')
    data = gt.Data()
    for frame in traj[::len(traj)//4000]:
        frame.wrap()
        data.add(frame)

    data.save(filename='GAP_rPBE0-D3_nvt_300K')
