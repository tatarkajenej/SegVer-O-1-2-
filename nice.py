#! /usr/bin/python3

from subprocess import run
from time import sleep

CERTIFICATE_DIR = "./cert/"
EXECUTABLE = "build/nice"
MAX_ATTEMPTS = 3
LOG_FILE = "./nice.log"
LINE_START = "The configuration is "
SUCCESS_WORD = "NON-DEFECTIVE"

def execute_task(which1, which2, u, v, filename=None):
    global log_handle

    log_handle.write("Starting task %i %i %i %i\n" % (which1, which2, u, v))
    if filename is None:
        filename = "%i-%i-%i-%i.run" % (which1, which2, u, v)
    filename = CERTIFICATE_DIR + filename

    # try several runs if you see defectivity (it might be a false negative due to the finite field arithmetic)
    success = False;
    attempts_done = 0;
    while not success and attempts_done < MAX_ATTEMPTS:
        run([EXECUTABLE, str(which1), str(which2), str(u), str(v), filename])
        with open(filename, "r") as cert:
            for line in cert:
                if len(line) > len(LINE_START) and line[:len(LINE_START)] == LINE_START:
                    success = (line[len(LINE_START):len(LINE_START)+len(SUCCESS_WORD)] == SUCCESS_WORD)
                    break
        attempts_done += 1
        if success:
            log_handle.write("-- attempt number %i succeeded\n" % attempts_done)
        else:
            log_handle.write("-- attempt number %i failed\n" % attempts_done)
            sleep(1)    # sleep one second to make sure the next attempt will have a different timestamp and thus a different random seed
    log_handle.write("\n")

configuration_names = {
    (0,0) : "A0",

    (1,0) : "B0",
    (1,1) : "B1",
    (1,2) : "B2",

    (2,0) : "C0",
    (2,1) : "C1",

    (3,0) : "D0",
    (3,1) : "D1",

    (4,0) : "E0",
    (4,1) : "E1",

    (5,0) : "F0",

    (6,0) : "G0"
}

tasks = [
    (6, 0, 12, 372),

    (5, 0, 11, 306),
    (5, 0, 10, 250),

    (4, 0, 9, 196),
    (4, 0, 9, 197),
    (4, 0, 8, 152),
    (4, 1, 8, 154),

    (3, 0, 7, 111),
    (3, 0, 7, 112),
    (3, 1, 7, 112),
    (3, 0, 6, 82),
    (3, 1, 6, 80),

    (2, 0, 5, 50),
    (2, 0, 5, 51),
    (2, 1, 5, 52),
    (2, 0, 4, 30),
    (2, 1, 4, 30),

    (1, 0, 3, 12),
    (1, 0, 3, 13),
    (1, 1, 3, 14),
    (1, 0, 2, 4),
    (1, 1, 2, 4),
    (1, 2, 2, 5),

    (0, 0, 3, 8),
    (0, 0, 3, 9),
    (0, 0, 4, 26)
]

# before starting, make sure the executable is up to date
run(["make", EXECUTABLE])

# execute sequentially
with open(LOG_FILE, "w") as log_handle:
    for args in tasks:
        execute_task(*args, "%s-%i-%i.run" % (configuration_names[args[:2]], args[2], args[3]))
