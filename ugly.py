#! /usr/bin/python3

from subprocess import run
from time import sleep

CERTIFICATE_DIR = "./cert/"
EXECUTABLE = "build/ugly"
MAX_ATTEMPTS = 3
LOG_FILE = "./ugly.log"
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
    (1,0) : "A1hat",
    (1,1) : "A1hat",
    (1,2) : "A2hat",
    (1,3) : "A2hat",
    (1,4) : "A3hat",
    (1,5) : "A3hat",

    (2,0) : "B1hat",
    (2,1) : "B1hat",
    (2,2) : "B2hat",
    (2,3) : "B2hat",
    (2,4) : "B3hat",
    (2,5) : "B3hat",

    (3,0) : "C1hat",
    (3,1) : "C1hat",
    (3,2) : "C2hat",
    (3,3) : "C2hat",

    (4,0) : "D1hat",
    (4,1) : "D1hat"
}

tasks = [
    (4, 0, 8, 21),
    (4, 1, 8, 23),

    (3, 0, 6, 15),
    (3, 1, 6, 17),
    (3, 2, 6, 19),
    (3, 3, 6, 21),

    (2, 0, 4, 9),
    (2, 1, 4, 11),
    (2, 2, 4, 13),
    (2, 3, 4, 15),
    (2, 4, 4, 17),
    (2, 5, 4, 19),

    (1, 0, 4, 9),
    (1, 1, 4, 11),
    (1, 2, 2, 7),
    (1, 3, 2, 9),
    (1, 4, 2, 11),
    (1, 5, 2, 13)
]

# before starting, make sure the executable is up to date
run(["make", EXECUTABLE])

# execute sequentially
with open(LOG_FILE, "w") as log_handle:
    for args in tasks:
        execute_task(*args, "%s-%i-%i.run" % (configuration_names[args[:2]], args[2], args[3]))
