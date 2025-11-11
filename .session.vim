let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/projects/trees
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
let s:shortmess_save = &shortmess
set shortmess+=aoO
badd +1 src/functions.py
badd +188 scripts/paper_plots.py
badd +42 config/config.yaml
badd +1 term://~/projects/trees//14279:/bin/bash
argglobal
%argdel
edit src/functions.py
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe '1resize ' . ((&lines * 34 + 36) / 73)
exe 'vert 1resize ' . ((&columns * 127 + 128) / 256)
exe '2resize ' . ((&lines * 35 + 36) / 73)
exe 'vert 2resize ' . ((&columns * 127 + 128) / 256)
exe 'vert 3resize ' . ((&columns * 128 + 128) / 256)
argglobal
balt config/config.yaml
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 22 - ((20 * winheight(0) + 17) / 34)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 22
normal! 040|
wincmd w
argglobal
if bufexists(fnamemodify("config/config.yaml", ":p")) | buffer config/config.yaml | else | edit config/config.yaml | endif
if &buftype ==# 'terminal'
  silent file config/config.yaml
endif
balt src/functions.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 42 - ((33 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 42
normal! 025|
wincmd w
argglobal
if bufexists(fnamemodify("term://~/projects/trees//14279:/bin/bash", ":p")) | buffer term://~/projects/trees//14279:/bin/bash | else | edit term://~/projects/trees//14279:/bin/bash | endif
if &buftype ==# 'terminal'
  silent file term://~/projects/trees//14279:/bin/bash
endif
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
let s:l = 12 - ((11 * winheight(0) + 35) / 70)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 12
normal! 02|
wincmd w
exe '1resize ' . ((&lines * 34 + 36) / 73)
exe 'vert 1resize ' . ((&columns * 127 + 128) / 256)
exe '2resize ' . ((&lines * 35 + 36) / 73)
exe 'vert 2resize ' . ((&columns * 127 + 128) / 256)
exe 'vert 3resize ' . ((&columns * 128 + 128) / 256)
tabnext 1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20
let &shortmess = s:shortmess_save
let &winminheight = s:save_winminheight
let &winminwidth = s:save_winminwidth
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
nohlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
