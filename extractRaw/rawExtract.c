#include <stdio.h>

#include "/Library/Frameworks/Mono.framework/Versions/Current/include/mono-2.0/mono/jit/jit.h"
#include "/Library/Frameworks/Mono.framework/Versions/Current/include/mono-2.0/mono/metadata/assembly.h"


int main (int argc, char *argv[]) {
    mono_config_parse (NULL);
    MonoDomain *domain;
    domain = mono_jit_init("extractRaw.exe");
    MonoAssembly *assembly;
    assembly = mono_domain_assembly_open (domain, "extractRaw.exe");
    if (!assembly) exit(1);
    //printf("Opened test_mono.exe\n");

    int retval;
    retval = mono_jit_exec (domain, assembly, argc, argv);
    mono_jit_cleanup (domain);

    return retval;
}
