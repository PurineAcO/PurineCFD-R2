"""Fluent 二进制 cas -> PurineCFD mesh.txt

规则来自 cases/oat15a/readrans.txt:
  * 段头形如 "(10 (<zone> <first> <last> ...)(" 之后的数字是 16 进制, 首尾编号含端点
  * 节点段 (10): 每记录 16 字节 = 2 个 double (x y); 若是 24 字节则取 3 个 double 的前两个
  * 面段   (13): 每记录 16 字节 = 4 个 int (n0 n1 c0 c1), c=0 表示无单元
  * 面类型由线程段 "(39 (<zone> <type> <name> <idx>)" 的 name 给出: interior / wall / far*
  * 输出 mesh.txt: 头 "N F C 3" + 三个组名 + (node)/(edge)/(cell) 三段
"""

import re
import struct
import sys

import numpy as np

CAS = sys.argv[1] if len(sys.argv) > 1 else 'cases/oat15a/1-21-00200.cas'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'cases/oat15a/mesh.txt'
INTER, WALL, FAR = 0, 1, 4


def sections(buf, tag, rec):
    """扫 '(<tag> (' 段, 用 (last-first+1)*rec 后是否紧跟 ')' 校验, 防止误命中数据体"""
    out = []
    pat = ('(' + str(tag) + ' (').encode()
    pos = 0
    while True:
        i = buf.find(pat, pos)
        if i < 0:
            return out
        pos = i + 1
        nl = buf.find(b'\n', i)
        if nl < 0:
            continue
        if buf[nl + 1:nl + 2] != b'(':
            continue
        toks = re.findall(rb'[0-9a-fA-F]+', buf[i + len(pat):nl])
        if len(toks) < 3:
            continue
        try:
            zone, first, last = (int(t, 16) for t in toks[:3])
        except ValueError:
            continue
        n = last - first + 1
        body = nl + 2
        if n <= 0 or body + n * rec >= len(buf):
            continue
        if buf[body + n * rec:body + n * rec + 1] != b')':
            continue
        out.append((zone, first, last, body, n))


def zone_kinds(buf):
    """(39 (<zone> <?> <name> <idx>) -> {zone: INTER/WALL/FAR}, 按名字里的关键字判定"""
    out = {}
    for m in re.finditer(rb'\(39 \(([^)]*)\)', buf):
        t = m.group(1).split()
        if len(t) < 2:
            continue
        try:
            zone = int(t[0])
        except ValueError:
            continue
        for tok in t[1:]:
            d = tok.decode(errors='ignore').strip().lower()
            if d == 'interior':
                out[zone] = INTER
                break
            if 'wall' in d:
                out[zone] = WALL
                break
            if 'far' in d:
                out[zone] = FAR
                break
    return out


def main():
    buf = open(CAS, 'rb').read()
    print('cas 大小 %.1f MB' % (len(buf) / 1e6))

    seg16 = sections(buf, 3010, 16)
    seg24 = sections(buf, 3010, 24)
    nseg = seg16 if sum(s[4] for s in seg16) >= sum(s[4] for s in seg24) else seg24
    rec = 16 if nseg is seg16 else 24
    if not nseg:
        print('Error: 未找到节点段 (10')
        return 1
    node = {}
    for zone, first, last, body, n in nseg:
        for k in range(n):
            x, y = struct.unpack_from('<2d', buf, body + rec * k)[:2]
            node[first + k] = (x, y)
    nid = {k: i + 1 for i, k in enumerate(sorted(node))}
    print('节点 %d (%d 段, 记录 %d 字节)' % (len(node), len(nseg), rec))

    kind = zone_kinds(buf)
    fseg = sections(buf, 3013, 16)
    if not fseg:
        print('Error: 未找到面段 (13')
        return 1
    face, zcount = [], {}
    for zone, first, last, body, n in fseg:
        t = kind.get(zone, INTER)
        zcount[zone] = zcount.get(zone, 0) + n
        for k in range(n):
            a, b, c0, c1 = struct.unpack_from('<4i', buf, body + 16 * k)
            face.append((nid.get(a, 0), nid.get(b, 0), c0, c1, t))
    print('面 %d (%d 段); 线程段: %s' % (len(face), len(fseg), kind))
    print('各 zone 面数: %s' % {z: zcount[z] for z in sorted(zcount)})

    ncell = max(max(f[2], f[3]) for f in face)
    cells = [[] for _ in range(ncell + 1)]
    for i, (a, b, c0, c1, t) in enumerate(face, 1):
        for c in (c0, c1):
            if c > 0:
                cells[c].append(i)
    bad = [c for c in range(1, ncell + 1) if len(cells[c]) != 4]
    print('单元 %d, 非四边单元 %d 个' % (ncell, len(bad)))
    if bad:
        print('  例: %s' % [(c, len(cells[c])) for c in bad[:5]])
    if any(f[0] == 0 or f[1] == 0 for f in face):
        print('警告: 有面引用了未在节点段出现的节点号')

    with open(OUT, 'w') as fp:
        fp.write('%d %d %d 3\n' % (len(node), len(face), ncell))
        for nm, kd in (('interior', 'INTER'), ('wall', 'WALL'), ('far', 'FAR')):
            fp.write('%s = %s\n' % (nm, kd))
        fp.write('(node)\n')
        for k in sorted(node):
            fp.write('%d %.16e %.16e\n' % (nid[k], node[k][0], node[k][1]))
        fp.write('(end)\n(edge)\n')
        for nm, tt in (('interior', INTER), ('wall', WALL), ('far', FAR)):
            fp.write('%s\n' % nm)
            cnt = 0
            for i, (a, b, c0, c1, t) in enumerate(face, 1):
                if t == tt:
                    fp.write('%d %d %d %d %d\n' % (i, a, b, c0, c1))
                    cnt += 1
            fp.write('(end)\n')
            print('组 %-9s %d 面' % (nm, cnt))
        fp.write('(end)\n(cell)\n')
        for c in range(1, ncell + 1):
            fp.write('%d %s\n' % (c, ' '.join(str(v) for v in sorted(cells[c]))))
    print('写出 %s' % OUT)

    wn = set()
    for a, b, c0, c1, t in face:
        if t == WALL:
            wn.add(a)
            wn.add(b)
    p = np.array([node[k] for k in sorted(node) if nid[k] in wn])
    if len(p) == 0:
        print('警告: 没有壁面节点, 无法给壁面几何')
        return 0
    span = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    i, j = np.unravel_index(span.argmax(), span.shape)
    print('--- OAT15A 壁面几何 ---')
    print('壁面节点 %d 个' % len(p))
    print('x 范围 %.6f .. %.6f  (跨度 %.6f)' % (p[:, 0].min(), p[:, 0].max(), p[:, 0].max() - p[:, 0].min()))
    print('y 范围 %.6f .. %.6f  (跨度 %.6f)' % (p[:, 1].min(), p[:, 1].max(), p[:, 1].max() - p[:, 1].min()))
    print('弦长(壁面最远两点距) = %.6f' % span[i, j])
    print('  LE = (%.6f, %.6f)   TE = (%.6f, %.6f)' % (p[i][0], p[i][1], p[j][0], p[j][1]))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
