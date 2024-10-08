
Reference:

https://git-scm.com/book/en/v2/Git-Branching-Remote-Branches

Bring data from remote server
```
git clone https://apatlpo@bitbucket.org/apatlpo/natl60_dimup.git
git log --oneline --decorate --graph --all
```

Create local branch from remote ones
```
git checkout -b ap_changes origin/ap_changes
git checkout -b sf_changes origin/sf_changes
```

Merge sf_changes into master
```
git checkout master
git merge sf_changes
git branch -d sf_changes
```

Now merge ap_changes into master
```
git merge ap_changes
(CONFLICT (content): Merge conflict in overview/plot_snapshot_2d_noproj.py)
vim overview/plot_snapshot_2d_noproj.py
git add overview/plot_snapshot_2d_noproj.py
git commit
```

Delete ap_changes branch
```
git branch -d ap_changes
```

Checkout a file from another branch
```
git checkout mybranch
git checkout otherbranch -- dev/file.py
```

