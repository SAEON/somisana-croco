The `*.py` function in this directory have been created as part of the somisana tools

The `croco_pytools` directory is included here as a git subtree, added as follows (run from the root dir of this repo):

```sh
git subtree add --prefix=crocotools_py/croco_pytools https://gitlab.inria.fr/croco-ocean/croco_pytools.git v1.0.1 --squash
```

Note we are presently cloning the v1.0.1 branch, which can be updated with further git subtree commands to update

Adding the `croco_pytools` in this way allows us to have access to all of their functions from within the somisana-croco repo, and also allows us to pull future updates made by the developers of this code 

