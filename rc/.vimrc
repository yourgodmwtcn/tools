set nocompatible
set modelines=0
set showmode
set nowrap
set cursorline

set tabstop=4
set expandtab
set gdefault
set autoindent

" Searching options
set ignorecase
set smartcase
set hlsearch
set incsearch

set title
set visualbell t_vb=
set noerrorbells
set wildchar=<TAB>
set wildmenu
set wildmode=list,longest
set ttyfast

set number
set hidden

set termencoding=utf-8
set encoding=utf-8
set lazyredraw                  " don't update the display while executing macros
set laststatus=2                " tell VIM to always put a status line in, even
                                "    if there is only one window
set cmdheight=2                 " use a status bar that is 2 rows high

nnoremap ; :
" Tells Vim to stop highlighting after ESC
nnoremap <CR> :noh<CR><CR>

set foldenable                  " enable folding
set foldcolumn=2                " add a fold column
set foldmethod=marker           " detect triple-{ style fold markers
set foldlevelstart=0            " start out with everything folded
set foldopen=block,hor,insert,jump,mark,percent,quickfix,search,tag,undo
                                " which commands trigger auto-unfold

" Treat ROMS input files as fortran files - syntax highlighting
au BufNewFile,BufRead *.in setfiletype fortran
au BufNewFile,BufRead *.h  setfiletype fortran

set bg=dark
set t_Co=256
colorscheme gardener

