import sys
sys.path.append(".")

from libcloudphxx import blk_1m

opts = blk_1m.opts_t()
opts.cond = True
opts.cevp = True
opts.revp = True
opts.conv = True
opts.accr = True
opts.sedi = False

rhod = [1]
th   = [270]
rv   = [0]
rc   = [0]
rr   = [0]
dt   = 1

blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
