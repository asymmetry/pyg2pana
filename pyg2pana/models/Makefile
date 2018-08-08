LIST := pbosted

all:
	@make -s -C $(LIST)

clean:
	@for item in $(LIST); do \
	    if [[ -d $$item ]]; then \
	        make clean -s -C $$item; \
	    fi; \
	done;