ifeq ($(STRUNA),1)

DEPS+=string
include $(ROOT)/string/plugin/linkobjs.inc

LOCALMACROS+=-D__STRUNA -include struna.defs
LOCALDEPS+=struna.defs

plugin_master.F90 : struna_vars.src struna_init.src struna_main.src struna_done.src
struna_vars.src : force_look
	sed -n '/!#vars 0/{:a;n;/!#end vars 0/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >$@
	sed -n '/!#vars 1/{:a;n;/!#end vars 1/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#vars 2/{:a;n;/!#end vars 2/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@

struna_init.src : force_look
	sed -n '/!#init 0/{:a;n;/!#end init 0/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >$@
	sed -n '/!#init 1/{:a;n;/!#end init 1/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#init 2/{:a;n;/!#end init 2/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#init 3/{:a;n;/!#end init 3/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	echo " if(allocated(sm_atomlist)) then ; do l=2,sm_atomlist(1)+1; ierr=atoms%uadd(sm_atomlist(l)) ; ierr=forced_atoms%uadd(sm_atomlist(l)) ; enddo ; endif ! make sure to add atoms to main ind arrays">>$@

struna_main.src : force_look
	sed -n '/!#main 0/{:a;n;/!#end main 0/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >$@
	sed -n '/!#main 1/{:a;n;/!#end main 1/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	echo " do i=2,sm_atomlist(1)+1; ierr=atoms%uadd(sm_atomlist(i)) ; ierr=forced_atoms%uadd(sm_atomlist(i)) ; enddo ; alist=master_compute_atomlist() ">>$@
	sed -n '/!#main 2/{:a;n;/!#end main 2/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#main 3/{:a;n;/!#end main 3/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#main 4/{:a;n;/!#end main 4/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@

struna_done.src : force_look
	sed -n '/!#done 0/{:a;n;/!#end done 0/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >$@
	sed -n '/!#done 1/{:a;n;/!#end done 1/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#done 2/{:a;n;/!#end done 2/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@
	sed -n '/!#done 3/{:a;n;/!#end done 3/b;p;ba}' $(ROOT)/string/plugin/sm_plugin.F90 >>$@

.SILENT : *src

endif
