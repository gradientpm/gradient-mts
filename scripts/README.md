# Rendering tools

This directory store all the scripts used for:
- generate scenes variations (different rendering alg.)
- launch over all these scene computations
- pack the results in a comprehensive way
- generate web page to understand these results

Moreover, files or directory named "igrida":
It is only for french cluster igrida. You can ignore these files.

## Setup

* Unpack the custom cygwin package
* Create a symbolic link to mitsuba.exe in /usr/local/bin. 
```
  ln -s <your mitsuba.exe>    /usr/local/bin/mitsuba
```
  This is for convenience. You can also use -m flag to point to Mitsuba executable. 

* Directory structure
```  
[home]
  gradient_pm_code
    scripts
      run.py 
      
  gradient_pm_data
    scenes
	references 
```
	
## Usage

* Generate XML files from the ori_*.xml template
```
python3 run.py -i ../../gradient_pm_data/scenes/ -o Test -f 0:5:0 -j 4 -l ./results/data/html_gpm.xml -d 5 -s cbox -g 
```

* Run Mitsuba 
```
python3 run.py -i ../../gradient_pm_data/scenes/ -o Test -f 0:5:0 -j 4 -l ./results/data/html_gpm.xml -d 5 -s cbox -m mitsuba -c
```

* Compact results
```
python3 run.py -i ../../gradient_pm_data/scenes/ -o Test -f 0:5:0 -j 4 -l ./results/data/html_gpm.xml -d 5 -s cbox -r ../../gradient_pm_data/references/ref_cbox.hdr -p
```

* Generate webpage
```
python3 run.py -i ../../gradient_pm_data/scenes/ -o Test -f 0:5:0 -j 4 -l ./results/data/html_gpm.xml -d 5 -s cbox -r ../../gradient_pm_data/references/ref_cbox.hdr -P
```


