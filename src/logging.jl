"""Logger that stores all messages in a `messages` field (array of strings).

The `in` operator can be used to check whether a string or regex is in any of
the recorded messages.

```julia
test_logger = QuantumTestLogger()
with_logger(test_logger) do
    #...
end
@test "message" in logger
```

This is a simplified version of `Logging.TestLogger`.
"""
struct QuantumTestLogger <: AbstractLogger
    messages::Vector{String}
    function QuantumTestLogger()
        new(String[])
    end
end

Logging.shouldlog(logger::QuantumTestLogger, level, _module, group, id) = true

Logging.min_enabled_level(logger::QuantumTestLogger) = Logging.Debug

Logging.catch_exceptions(logger::QuantumTestLogger) = false

function Logging.handle_message(
    logger::QuantumTestLogger,
    level::Logging.LogLevel,
    message,
    args...;
    kwargs...
)
    levelstr = string(level)
    push!(logger.messages, "$levelstr: $message")
    nothing
end

function Base.in(s, logger::QuantumTestLogger)
    return any([occursin(s, msg) for msg in logger.messages])
end
