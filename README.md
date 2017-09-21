
## Install

Windows:

```powershell
cd <compile path>
cmake -G "Visual Studio 10 2010 Win64" <path/to/neighbor_based> -DINSTALL_PREFIX=<path/to/bin> -DARCH=64
```



## subtree 更新

```shell
# pull
git fetch utilsclass master
git subtree pull --prefix=src/UtilsClass utilsclass master --squash
# push
git subtree push --prefix=src/UtilsClass utilsclass master
```





