#
# Makefile for mmf
#

include ../makelist

SRCS = 	mmf.c parse_args.c alloc_space.c build_lists.c \
	setup_cont.c decl_control.c control_addr.c \
	control_var.c read_params.c sort_dims.c sort_params.c sort_vars.c \
	var_addr.c declvar.c str_to_vals.c\
	declparam.c param_addr.c getdim.c timing.c getparam.c umalloc_etc.c \
	julday.c getvar.c julconvert.c readvar.c decldim.c \
	get_times.c batch_run.c read_control.c dim_addr.c reset_dim.c read_line.c \
	get_elem_add.c read_vars.c getdimname.c \
	save_params.c load_param.c check_vars.c \
	create_vstats.c free_vstats.c write_vstats.c \
	call_modules.c call_setdims.c \
	read_datainfo.c putvar.c print_params.c print_vars.c \
	print_model_info.c batch_run_functions.c graph_single_run.c \
	control_array.c call_setdims.c call_modules.c


MMFOBJS = ${SRCS:.c=.o}

.c.o:
	$(CC) $(CFLAGS) -c $<

#
# Standard Targets for Users
#
all: $(MMFLIB)

$(MMFLIB): $(MMFOBJS)
# Create lib directory, if necessary
	@if [ ! -d $(LIBDIR) ]   ; then        \
   	   mkdir $(LIBDIR) ;                   \
	   echo  Created directory $(LIBDIR) ; \
	fi
	$(AR) $(MMFLIB) $(MMFOBJS)
	$(RANLIB) $(MMFLIB)

clean:
	$(RM) $(MMFLIB) $(MMFOBJS) *~
