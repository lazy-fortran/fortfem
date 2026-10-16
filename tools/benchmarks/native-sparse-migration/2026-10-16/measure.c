#define _GNU_SOURCE
#include <sched.h>
#include <spawn.h>
extern char **environ;
#include <sys/resource.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

static double seconds(struct timeval t) { return t.tv_sec + 1e-6*t.tv_usec; }
int main(int argc, char **argv) {
    if (argc != 3) return 2;
    int log = open(argv[2], O_CREAT|O_WRONLY|O_TRUNC, 0600);
    if (log < 0) return 3;
    cpu_set_t cpus; CPU_ZERO(&cpus); CPU_SET(0, &cpus);
    if (sched_setaffinity(0, sizeof(cpus), &cpus) != 0) return 6;
    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_adddup2(&actions, log, STDOUT_FILENO);
    posix_spawn_file_actions_adddup2(&actions, log, STDERR_FILENO);
    posix_spawn_file_actions_addclose(&actions, log);
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);
    pid_t child;
    char *child_argv[] = {argv[1], NULL};
    if (posix_spawn(&child, argv[1], &actions, NULL, child_argv, environ) != 0)
        return 4;
    posix_spawn_file_actions_destroy(&actions);
    close(log);
    int status; struct rusage usage;
    while (wait4(child, &status, 0, &usage) < 0) if (errno != EINTR) return 5;
    clock_gettime(CLOCK_MONOTONIC, &end);
    int exit_status = WIFEXITED(status) ? WEXITSTATUS(status) : 128+WTERMSIG(status);
    printf("{\"wall_s\":%.9f,\"user_s\":%.9f,\"system_s\":%.9f,"
           "\"cpu_s\":%.9f,\"maxrss_kib\":%ld,\"exit_status\":%d}\n",
           end.tv_sec-begin.tv_sec + 1e-9*(end.tv_nsec-begin.tv_nsec),
           seconds(usage.ru_utime), seconds(usage.ru_stime),
           seconds(usage.ru_utime)+seconds(usage.ru_stime), usage.ru_maxrss, exit_status);
    return exit_status;
}
