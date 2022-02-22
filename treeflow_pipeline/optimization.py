import tensorflow as tf


class RobustOptimizer:  # Pseudo-optimizer
    def __init__(self, inner, max_retries=10):  # TODO: Count number of failed steps
        self.inner = inner
        self.max_retries = max_retries
        self.retries = tf.Variable(0)

    def apply_gradients(self, grads_and_vars, name=None, **kwargs):
        """Apply gradients to variables.

        Args:
        grads_and_vars: List of (gradient, variable) pairs as returned by
            `compute_gradients()`.
        global_step: Optional `Variable` to increment by one after the
            variables have been updated.
        name: Optional name for the returned operation.  Default to the
            name passed to the `Optimizer` constructor.
        """
        grads_and_vars = list(grads_and_vars)
        grads = [grad for grad, var in grads_and_vars]
        # all_finite = tf.reduce_all([tf.reduce_all(tf.math.is_finite(x)) for x in grads])
        # return tf.cond(
        #     all_finite,
        #     lambda: self.inner.apply_gradients(grads_and_vars, name=name, **kwargs),
        #     lambda: tf.constant(True, name=name),
        # )
        any_nan = tf.reduce_any([tf.reduce_any(tf.math.is_nan(x)) for x in grads])

        def nan_handler():
            assertion = tf.assert_less(self.retries, self.max_retries)
            with tf.control_dependencies([assertion]):
                self.retries.assign_add(1)
            return tf.no_op()

        def update_handler():
            self.retries.assign(0)
            return self.inner.apply_gradients(grads_and_vars, name=name, **kwargs)

        return tf.cond(
            any_nan,
            nan_handler,
            update_handler,
        )