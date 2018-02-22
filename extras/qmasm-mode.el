;;; qmasm-mode.el --- Emacs mode for editing QMASM code

;; Author: Scott Pakin <pakin@lanl.gov>
;; Keywords: tools, languages

(defconst qmasm-mode-syntax-table
  (let ((table (make-syntax-table)))
    ;; '"' is a string delimiter.
    (modify-syntax-entry ?\" "\"" table)

    ;; Comments go from "#" to the end of the line.
    (modify-syntax-entry ?\# "<" table)
    (modify-syntax-entry ?\n ">" table)
    table))

(defvar qmasm-highlights
  ;; Directives are treated as functions.
  '(("![a-z_.]+" . font-lock-function-name-face)))

(define-derived-mode qmasm-mode prog-mode "QMASM"
  :syntax-table qmasm-mode-syntax-table
  (setq font-lock-defaults '(qmasm-highlights))
  (font-lock-fontify-buffer))
