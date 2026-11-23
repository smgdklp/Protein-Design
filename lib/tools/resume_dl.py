#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
可续传下载器（针对只支持 open-ended Range 的站点，例如 dln1.ncdn.ec）

用法：
    python resume_dl.py "<下载链接>" "D:\\some\\out.pdf"

特性：
  * 断点续传：复用同名 .part，本次用 Range: bytes=<已下载字节>- 接着下
  * 只支持 open-ended Range 的站点也能用（该站闭区间会返回 416）
  * 链接过期（410/403）时立刻停下并保留 .part，换新链接后可继续
  * 每 2 秒打印进度/速度/预计剩余时间
  * 结束后校验：字节数、%PDF 头、%%EOF 尾
"""
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.request

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/121 Safari/537.36"
CHUNK = 1 << 16


def ctx():
    c = ssl.create_default_context()
    c.check_hostname = False
    c.verify_mode = ssl.CERT_NONE
    return c


def total_size(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Range": "bytes=0-"})
    with urllib.request.urlopen(req, timeout=45, context=ctx()) as r:
        cr = r.headers.get("Content-Range")
        if cr:
            m = re.search(r"/(\d+)$", cr)
            if m:
                return int(m.group(1))
        cl = r.headers.get("Content-Length")
        return int(cl) if cl else None


def head_tail_ok(path):
    with open(path, "rb") as f:
        head = f.read(5)
        f.seek(max(0, os.path.getsize(path) - 64))
        tail = f.read()
    return head[:4] == b"%PDF", b"%%EOF" in tail


def download(url, out, max_retry=999, stall_retry=6):
    part = out + ".part"
    have = os.path.getsize(part) if os.path.exists(part) else 0
    try:
        total = total_size(url)
    except urllib.error.HTTPError as e:
        if e.code in (403, 410):
            print("链接已过期或被拒（HTTP %d）。请重新获取下载链接后再次运行，"
                  "本次已下载的 %d 字节会保留在 .part 里。" % (e.code, have))
            return 1
        raise
    if total is None:
        print("拿不到总大小，无法续传。")
        return 1
    print("总大小: %.1f MB，已有: %.1f MB，待下载: %.1f MB"
          % (total / 1048576, have / 1048576, (total - have) / 1048576))
    if have >= total:
        print("已完整。")
        return finish(part, out)

    retry = 0
    while have < total and retry < max_retry:
        req = urllib.request.Request(url, headers={
            "User-Agent": UA,
            "Range": "bytes=%d-" % have,
        })
        t0 = time.time()
        got = 0
        try:
            with urllib.request.urlopen(req, timeout=45, context=ctx()) as r:
                if r.status == 200 and have > 0:
                    print("服务器忽略了 Range（返回 200），只能从头下。把 .part 删了重来。")
                    return 1
                with open(part, "ab") as f:
                    last = time.time()
                    while True:
                        b = r.read(CHUNK)
                        if not b:
                            break
                        f.write(b)
                        got += len(b)
                        have += len(b)
                        if time.time() - last >= 2:
                            last = time.time()
                            sp = got / max(1e-6, time.time() - t0) / 1024
                            left = (total - have) / max(1e-6, sp * 1024)
                            print("  %.1f%%  %6.1f KB/s  剩余约 %s"
                                  % (have * 100.0 / total, sp, fmt(left)), flush=True)
            retry = 0
        except urllib.error.HTTPError as e:
            if e.code in (403, 410):
                print("链接过期/被拒（HTTP %d）。已下 %d 字节保留，换新链接后重跑本命令即可续传。"
                      % (e.code, have))
                return 1
            retry += 1
            print("HTTP %d，2s 后重试（第 %d 次）" % (e.code, retry))
            time.sleep(2)
        except Exception as e:
            retry += 1
            print("中断（%s），2s 后重试（第 %d 次）" % (type(e).__name__, retry))
            time.sleep(2)
        if got == 0:
            stall_retry -= 1
            if stall_retry <= 0:
                print("连续多次没进展，可能链路被限速。建议开代理/加速器后重跑（会续传）。")
                return 1
    return finish(part, out)


def fmt(s):
    s = int(s)
    return "%d:%02d:%02d" % (s // 3600, s % 3600 // 60, s % 60)


def finish(part, out):
    ok_head, ok_tail = head_tail_ok(part)
    print("校验: %%PDF=%s  %%%%EOF=%s" % (ok_head, ok_tail))
    if os.path.exists(out):
        os.remove(out)
    os.replace(part, out)
    print("完成:", out, "%.1f MB" % (os.path.getsize(out) / 1048576))
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(download(sys.argv[1], sys.argv[2]))
