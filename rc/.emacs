(setq column-number-mode t)

;;more real estate
(tool-bar-mode -1)
(menu-bar-mode -1)

;;line number column
(global-linum-mode)
(setq linum-format "%d ")

(add-to-list 'load-path "~/.emacs.d/matlab-emacs/matlab-emacs/")
(load-library "matlab-load")

(global-font-lock-mode t)
(setq font-lock-maximum-decoration t)