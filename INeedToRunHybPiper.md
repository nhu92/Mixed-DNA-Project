# Wish I could rerun the HybPiper!
This is a note-style markdown for documenting steps for a rerun of HybPiper on Mixed DNA samples

1. Delete conda, re-install conda
I had to delete conda in the HPCC of TTU since it became super slow for the past couple of years. Now I am graduated and temporarily free from worrying about reproducing my old analysis (I hope) (And I am working on documenting those pipelines, slowly. See [HomologousAnnotationPipeline](https://github.com/gudusanjiao/HomologousAnnotationPipeline)). As suggested on Anaconda website, there are two ways to uninstall conda from remote computer clusters - the full uninstall and the 'nearly' full uninstall. The full uninstall is automatically operated by conda but it requires a functional conda (which I don't have and that's why I want to reinstall right?). Thus, I did a manual uninstall, which is basically just deleting everything from conda directories and wiping out most of the traces from user files.
As suggested, here are command lines and additional steps to take into action.
```bash
rm -rf ~/anaconda3
rm -rf ~/conda
```
Also, search for any "suspicious" files created by conda in `/home/<user>`.

After this, search for files that look like `.bash*` (for example: `.bashrc`). It requires to use `ls -a` to show those hidden files. Delete anything that popped up in those files that you think relates to conda.

Now, I get a fresh and clean "conda-free" environment for my HPCC account. It's about time to restart everything!


