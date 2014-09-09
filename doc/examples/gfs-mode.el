;;; gfs-mode.el

;;; Copyright: (C) 2010 Stephane Popinet
;; 
;;     This program is free software; you can redistribute it and/or
;;     modify it under the terms of the GNU General Public License as
;;     published by the Free Software Foundation; either version 2 of
;;     the License, or (at your option) any later version.
;;     
;;     This program is distributed in the hope that it will be useful,
;;     but WITHOUT ANY WARRANTY; without even the implied warranty of
;;     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
;;     GNU General Public License for more details.
;;     
;;     You should have received a copy of the GNU General Public License
;;     along with GNU Emacs; if not, write to the Free Software
;;     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;     02110-1301 USA
;;
;; To use this package, you can save this file somewhere in your
;; load-path and put the following in your .emacs at a minimum:
;;
;; (require 'gfs-mode)

(define-derived-mode gfs-mode shell-script-mode "Gerris"
  "Major mode for editing Gerris simulation files."
  
  (require 'gfs-keywords)

  (defvar gfs-browse-base "http://gfs.sourceforge.net/wiki/index.php/"
    "First part of URL used to display documentation on the Gerris website.")

  (defvar gfs-ref-regexp
    (eval-when-compile
      (concat "\\<" (regexp-opt gfs-abbrevs t) "\\>"))
    "Regular expression compiled Gerris keywords.")

  (defvar gfs-modules-regexp
    (eval-when-compile
      (concat "\\<" (regexp-opt gfs-modules t) "\\>"))
    "Regular expression compiled Gerris modules.")
  
  (define-key gfs-mode-map [mouse-2] 'gfs-mode-mouse-2)
  (define-key gfs-mode-map [follow-link] 'mouse-face)

  (defun gfs-clickable-refs (limit)
    "Font-lock function which finds Gerris keywords and makes them clickable."
    (if	(re-search-forward (eval gfs-ref-regexp) limit t)
	(progn
	  (add-text-properties (match-beginning 0) (match-end 0)
			       (list 'mouse-face 'highlight
				     'gfs-keyword (match-string 0)
				     'help-echo "mouse-2: documentation"
				     'rear-nonsticky '(mouse-face gfs-keyword help-echo)))
	  t)))

  (defun gfs-clickable-modules (limit)
    "Font-lock function which finds Gerris modules and makes them clickable."
    (if	(re-search-forward (eval gfs-modules-regexp) limit t)
	(progn
	  (add-text-properties (match-beginning 0) (match-end 0)
			       (list 'mouse-face 'highlight
				     'gfs-module (match-string 0)
				     'help-echo "mouse-2: documentation"
				     'rear-nonsticky '(mouse-face gfs-module help-echo)))
	  t)))

  (defun gfs-comments (limit)
    "Font-lock function which finds Gerris comments."
    (re-search-forward "#.*$" limit t))

  (defconst gfs-font-lock-keywords
    (list 
     '(gfs-clickable-refs (0 'font-lock-function-name-face t))
     '(gfs-clickable-modules (0 'font-lock-type-face t))
     '(gfs-comments (0 'font-lock-comment-face t)))
    "Font-lock-keywords to be added when gfs-mode is active.")

  (defun gfs-url-create (ref-string module)
    "Returns REF-STRING without carriage returns and with spaces converted
to + signs, useful when creating a URL to lookup on the Gerris website."
    (with-temp-buffer
      (insert gfs-browse-base)
      (if module
	  (progn 
	    (insert "Object_hierarchy#")
	    (insert (capitalize ref-string)))
	(progn 
	  (unless (string= (substring ref-string 0 3) "Gfs")
	    (insert "Gfs"))
	  (insert ref-string)))
      (buffer-string)))

  (defun gfs-browse-reference (reference &optional module)
    "Wrapper function to call standard Emacs browser function for REFERENCE."
    (message "Linking to Gerris website for %s..." reference)
    (browse-url (gfs-url-create reference module)))
  
  (defun gfs-mode-mouse-2 (event arg)
    "Fetch documentation for keyword under the mouse click."
    (interactive "e\nP")
    (let (my-keyword)
      (save-excursion
	(set-buffer (window-buffer (posn-window (event-end event))))
	(goto-char (posn-point (event-end event)))
	(setq my-keyword (get-text-property (point) 'gfs-keyword)))
      (if my-keyword
	  (progn
	    (select-window (posn-window (event-end event)))
	    (gfs-browse-reference my-keyword))
	(progn
	  (setq my-keyword (get-text-property (point) 'gfs-module))
	  (if my-keyword
	      (progn
		(select-window (posn-window (event-end event)))
		(gfs-browse-reference my-keyword t))
	    (mouse-yank-at-click event arg)
	    )))))

  (font-lock-add-keywords nil gfs-font-lock-keywords)

  ;; load keywords for autocompletion with dabbrev
  (find-file-noselect (locate-file "gfs-keywords.el" load-path) t)
  (setq case-fold-search nil)

  (column-number-mode 1)
)

(add-to-list 'auto-mode-alist '("\\.gfs\\'" . gfs-mode))

(provide 'gfs-mode)
