grep -o -n --colour=never '\\beginobject{.*}' ../src/*.[ch] | awk 'BEGIN{FS="{|}|:"}{print $1 ":" $2,$4}'
