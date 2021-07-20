#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vi:ts=4:et
#
#
# Usage: python retriever-multi.py <file with URLs to fetch> [<# of
#          concurrent connections>]
#
# Modified version of: https://github.com/pycurl/pycurl/blob/master/examples/retriever-multi.py
#

import sys
import os
import pycurl
from _include import log

# We should ignore SIGPIPE when using pycurl.NOSIGNAL - see
# the libcurl tutorial for more info.
try:
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_IGN)
except ImportError:
    pass


# Get args
num_conn = 10
try:
    if sys.argv[1] == "-":
        urls = sys.stdin.readlines()
    else:
        urls = open(sys.argv[1]).readlines()
    if len(sys.argv) == 3:
        num_conn = int(sys.argv[2])
except:
    log("Usage: %s <file with URLs to fetch> [<# of concurrent connections>]" % sys.argv[0])
    raise SystemExit


# Make a queue with (url, filename) tuples
queue = []
for url in urls:
    url = url.strip()
    if not url or url[0] == "#":
        continue
    # filename = "doc_%03d.dat" % (len(queue) + 1)
    filename = os.path.basename(url)
    queue.append((url, filename))


# Check args
assert queue, "no URLs given"
num_urls = len(queue)
num_conn = min(num_conn, num_urls)
assert 1 <= num_conn <= 10000, "invalid number of concurrent connections"
log("PycURL %s (compiled against 0x%x)" % (pycurl.version, pycurl.COMPILE_LIBCURL_VERSION_NUM))
log("----- Getting", num_urls, "URLs using", num_conn, "connections -----")


# Pre-allocate a list of curl objects
m = pycurl.CurlMulti()
m.handles = []
for i in range(num_conn):
    c = pycurl.Curl()
    c.fp = None
    c.setopt(pycurl.FOLLOWLOCATION, 1)
    c.setopt(pycurl.MAXREDIRS, 5)
    c.setopt(pycurl.CONNECTTIMEOUT, 30)
    c.setopt(pycurl.TIMEOUT, 300)
    c.setopt(pycurl.NOSIGNAL, 1)
    m.handles.append(c)


# Main loop
freelist = m.handles[:]
num_processed = 0
num_failed = 0
while num_processed < num_urls:
    # Check if the file aready exist, and if it does check the file size (to do)
    if num_processed % 10 == 0 or num_processed == num_urls-1:
        log('\r- processed % 5d out of % 5d ( % 4d failed)' % (num_processed, num_urls, num_failed), end='')
    if len(queue) < 1:
        break
    # If file exists, don't download it again. THIS DOES NOT CHECK INCOMPLETE DOWNLOADS!
    if os.path.isfile(queue[0][1]):
        num_processed += 1
        del queue[0]
        continue
    # If there is an url to process and a free curl object, add to multi stack
    while queue and freelist:
        url, filename = queue.pop(0)
        c = freelist.pop()
        c.fp = open(filename, "wb")
        c.setopt(pycurl.URL, url)
        c.setopt(pycurl.WRITEDATA, c.fp)
        m.add_handle(c)
        # store some info
        c.filename = filename
        c.url = url
    # Run the internal curl state machine for the multi stack
    while 1:
        ret, num_handles = m.perform()
        if ret != pycurl.E_CALL_MULTI_PERFORM:
            break
    # Check for curl objects which have terminated, and add them to the freelist
    while 1:
        num_q, ok_list, err_list = m.info_read()
        for c in ok_list:
            c.fp.close()
            c.fp = None
            m.remove_handle(c)
            # print("Success:", c.filename, c.url, c.getinfo(pycurl.EFFECTIVE_URL))
            freelist.append(c)
        for c, errno, errmsg in err_list:
            c.fp.close()
            c.fp = None
            m.remove_handle(c)
            # print("Failed: ", c.filename, c.url, errno, errmsg)
            freelist.append(c)
        num_processed = num_processed + len(ok_list) + len(err_list)
        num_failed = num_failed + len(err_list)
        if num_q == 0:
            break
    # Currently no more I/O is pending, could do something in the meantime
    # (display a progress bar, etc.).
    # We just call select() to sleep until some more data is available.
    m.select(1.0)


# Cleanup
for c in m.handles:
    if c.fp is not None:
        c.fp.close()
        c.fp = None
    c.close()
m.close()
log('\r- processed % 5d out of % 5d ( % 4d failed)' % (num_processed, num_urls, num_failed))
