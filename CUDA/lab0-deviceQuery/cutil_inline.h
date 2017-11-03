#define cutilExit(argc, argv) __cutilExit(argc, argv)

inline void __cutilExit(int argc, char **argv) {
    if (!checkCmdLineFlag(argc, (const char**)argv, "noprompt")) {
        printf("\nPress ENTER to exit...\n");
        fflush( stdout);
        fflush( stderr);
        getchar();
    }
    cudaDeviceReset();
    exit(EXIT_SUCCESS);
}

